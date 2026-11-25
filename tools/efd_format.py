"""Small, non-evaluating parser/compiler for EFD arithmetic expressions."""
import ast
import re
import keyword
from functools import lru_cache

TOKEN = re.compile(r'\s*(\*\*|[A-Za-z_][A-Za-z_0-9]*|[0-9]+|[()+*/^=-])')


def explicit(expression):
    tokens = []; pos = 0
    while pos < len(expression):
        match = TOKEN.match(expression, pos)
        if not match:
            if expression[pos:].strip() == '': break
            raise ValueError('Unsupported arithmetic: ' + expression)
        tokens.append(match.group(1)); pos = match.end()
    result = []
    for token in tokens:
        if result and (result[-1] == ')' or re.fullmatch(r'[A-Za-z_0-9]+', result[-1])) and (
                token == '(' or re.fullmatch(r'[A-Za-z_0-9]+', token)):
            result.append('*')
        result.append(token)
    return ''.join(result).replace('**', '^')


def normalize(text):
    aliases={name:'p_'+name for name in keyword.kwlist if re.search(r'\b'+name+r'\b',text)}
    for line in text.splitlines():
        if line.startswith('parameter '): name=line.split()[1]
        elif line.startswith('compute '): name=line[8:].split('=',1)[0].strip()
        else: continue
        if re.fullmatch(r'[0-9]+[A-Za-z_][A-Za-z_0-9]*',name): aliases[name]='p_'+name
    result=[]
    for line in text.splitlines():
        tag, sep, value = line.partition(' ')
        if tag not in ('source','name'):
            for old,new in aliases.items():value=re.sub(r'\b'+re.escape(old)+r'\b',new,value)
            line=tag+sep+value
        if tag in {'satisfying','addition','doubling','tripling','negation','neutral',
                   'toweierstrass','fromweierstrass','assume','compute','toaffine','tosystem','homogweight'}:
            if '=' not in value: raise ValueError(line)
            lhs,rhs=value.split('=',1)
            equal=' == ' if tag=='satisfying' and line.startswith('satisfying ') and '==' in line else ' = '
            if rhs.startswith('='): rhs=rhs[1:]
            if tag=='satisfying' and lhs.strip() in {'x^3+y^3+1','y^2','a*x^3+y^3+1'}: equal=' == '
            line=tag+' '+explicit(lhs.strip())+equal+explicit(rhs.strip())
        elif re.match(r'a[013246] =',line):
            lhs,rhs=line.split('=',1);line=lhs.strip()+' = '+explicit(rhs)
        result.append(line.rstrip())
    return '\n'.join(result).rstrip()+'\n'


def compile_op3(computations):
    """Lower arithmetic ASTs to binary/unary operations, without eval or simplification."""
    code=[]; counter=0
    used=set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*','\n'.join(computations)))
    def temp():
        nonlocal counter
        while True:
            name='efdtmp'+str(counter);counter+=1
            if name not in used:used.add(name);return name
    def lower(node, target=None):
        if isinstance(node,ast.Name): value=node.id
        elif isinstance(node,ast.Constant) and type(node.value)==int: value=str(node.value)
        elif isinstance(node,ast.BinOp) and type(node.op) in {ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow}:
            left=lower(node.left);right=lower(node.right)
            value=left+{ast.Add:'+',ast.Sub:'-',ast.Mult:'*',ast.Div:'/',ast.Pow:'^'}[type(node.op)]+right
        elif isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub): value='-'+lower(node.operand)
        else: raise ValueError('Unsupported AST '+ast.dump(node))
        if target is not None: code.append(target+' = '+value);return target
        if isinstance(node,(ast.Name,ast.Constant)):return value
        name=temp();code.append(name+' = '+value);return name
    for assignment in computations:
        node=ast.parse(assignment.replace('^','**'),mode='exec').body
        if len(node)!=1 or not isinstance(node[0],ast.Assign) or len(node[0].targets)!=1 or not isinstance(node[0].targets[0],ast.Name):raise ValueError(assignment)
        lower(node[0].value,node[0].targets[0].id)
    return '\n'.join(code)+'\n'


@lru_cache(maxsize=32768)
def expression_tree(expression):
    return ast.parse(expression.replace('^','**'),mode='eval').body


def evaluate(expression, env, p):
    """Evaluate a trusted-format AST using exact finite-field operations, never Python eval."""
    def run(node):
        if isinstance(node,ast.Name):return env[node.id] % p
        if isinstance(node,ast.Constant) and type(node.value)==int:return node.value % p
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -run(node.operand)%p
        if isinstance(node,ast.BinOp):
            u=run(node.left)
            if isinstance(node.op,ast.Pow):
                if not isinstance(node.right,ast.Constant) or type(node.right.value)!=int:raise ValueError('Nonconstant exponent')
                return pow(u,node.right.value,p)
            v=run(node.right)
            if isinstance(node.op,ast.Add):return (u+v)%p
            if isinstance(node.op,ast.Sub):return (u-v)%p
            if isinstance(node.op,ast.Mult):return u*v%p
            if isinstance(node.op,ast.Div):return u*pow(v,-1,p)%p
        raise ValueError(ast.dump(node))
    return run(expression_tree(expression))
