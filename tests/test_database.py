"""Integrity and arithmetic-domain regression tests for the expanded database."""
import ast
import unittest
from tools.efd_format import normalize,compile_op3,evaluate
from tools.verify_binary_imports import GF2,ev
from tools.validate_database import validate


class Database(unittest.TestCase):
    def test_index_and_source_programs(self):
        counts=validate()
        self.assertGreaterEqual(counts['formula_entries'],639)
        self.assertEqual(counts['binary_model_families'],3)

    def test_integer_constants_and_binary_field_elements_are_distinct(self):
        F=GF2(19)
        self.assertEqual(ev('2',{},F),0)
        self.assertEqual(ev('a',{'a':2},F),2)
        self.assertEqual(F.mul(8,2),3)  # x^4 reduces to x+1 modulo x^4+x+1
        self.assertEqual(ev('a+a',{'a':7},F),0)
        for modulus in (19,37,67):
            field=GF2(modulus)
            for a in range(1,field.q):self.assertEqual(field.mul(a,field.inv(a)),1)
        self.assertEqual(evaluate('a+a',{'a':7},101),14)

    def test_parser_preserves_numeric_names_and_keyword_names(self):
        text=normalize('parameter 3overd\nassume 3overd*d = 3\ncompute 2P = 2 P\ncompute lambda = (2P+1)*3overd\n')
        self.assertIn('parameter p_3overd',text)
        self.assertIn('compute p_2P = 2*P',text)
        code=compile_op3([l[8:] for l in text.splitlines() if l.startswith('compute ')])
        ast.parse(code.replace('^','**'))

    def test_ast_evaluator_rejects_non_arithmetic(self):
        for expression in ('__import__("os")','a.real','[a][0]'):
            with self.assertRaises((ValueError,SyntaxError,KeyError)):evaluate(expression,{'a':1},101)


if __name__=='__main__':unittest.main()
