import ast
import operator as op
import yaml
from pathlib import Path

_ALLOWED_BINOPS = {
    ast.Add: op.add,
    ast.Sub: op.sub,
    ast.Mult: op.mul,
    ast.Div: op.truediv,
    ast.Pow: op.pow,
}

_ALLOWED_UNARYOPS = {
    ast.UAdd: op.pos,
    ast.USub: op.neg,
}


class UserExpression:
    def __init__(self, expression: str):
        self.expression = expression
        self._tree = ast.parse(expression, mode="eval")

    def value(self, context: dict) -> float:
        return float(self._eval_node(self._tree.body, context))

    def _eval_node(self, node, context):
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (int, float)):
                return node.value
            raise TypeError(f"Unsupported constant: {node.value!r}")

        if isinstance(node, ast.BinOp):
            op_type = type(node.op)
            if op_type not in _ALLOWED_BINOPS:
                raise TypeError(f"Operator not allowed: {op_type.__name__}")
            left = self._eval_node(node.left, context)
            right = self._eval_node(node.right, context)
            return _ALLOWED_BINOPS[op_type](left, right)

        if isinstance(node, ast.UnaryOp):
            op_type = type(node.op)
            if op_type not in _ALLOWED_UNARYOPS:
                raise TypeError(f"Unary operator not allowed: {op_type.__name__}")
            operand = self._eval_node(node.operand, context)
            return _ALLOWED_UNARYOPS[op_type](operand)

        if isinstance(node, ast.Name):
            return self._resolve_name(node.id, context)

        if isinstance(node, ast.Attribute):
            return self._resolve_attribute_chain(node, context)

        raise TypeError(f"Unsupported syntax: {ast.dump(node)}")

    def _resolve_name(self, name: str, context: dict):
        if name not in context:
            raise KeyError(f"Unknown name: {name}")
        return context[name]

    def _resolve_attribute_chain(self, node: ast.Attribute, context: dict):
        parts = []
        current = node

        while isinstance(current, ast.Attribute):
            parts.append(current.attr)
            current = current.value

        if not isinstance(current, ast.Name):
            raise TypeError("Only simple dotted paths are allowed")

        parts.append(current.id)
        parts.reverse()

        value = context
        for part in parts:
            if not isinstance(value, dict):
                raise KeyError(f"Cannot descend into non-dict at: {part}")
            if part not in value:
                raise KeyError(f"Unknown path: {'.'.join(parts)}")
            value = value[part]

        return value
