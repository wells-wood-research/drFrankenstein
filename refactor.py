import ast
import os
import re

# --- Naming Convention Helpers ---
def to_snake_case(name):
    if not name or name == "_" or (name.startswith("__") and name.endswith("__")):
        return name
    if name.isupper() and len(name) > 1: # If all upper (e.g. URL, ID), convert to lower
        return name.lower()
    # Handle cases like "CAMELCase" -> "camel_case"
    s0 = re.sub('([A-Z]+)([A-Z][a-z])', r'\1_\2', name) # Handle sequences of caps followed by cap-lower: IPAddress -> IP_Address
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', s0)    # Basic camel to snake: myVar -> my_Var
    s2 = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower() # Lowercase and final split: my_Var -> my_var, IP_Address -> ip_address
    return s2.replace("__", "_")

def to_camel_case(name):
    if not name or name == "_" or (name.startswith("__") and name.endswith("__")):
        return name

    leading_underscores = ""
    idx = 0
    while idx < len(name) and name[idx] == '_':
        leading_underscores += "_"
        idx += 1
    name = name.lstrip('_')
    if not name: return leading_underscores # Was all underscores

    if '_' not in name: # No underscores
        if name[0].islower() and (len(name) == 1 or any(c.isupper() for c in name[1:]) or not any(c.isupper() for c in name)): # Already camelCase e.g. myVar or alllower
            return leading_underscores + name
        if name[0].isupper(): # PascalCase e.g. MyVar or ALLCAPS e.g. URL
            # Convert MyVar to myVar, URL to url, MyURL to myUrl
            # Find first sequence of uppercase letters that is followed by a lowercase, or is the whole string
            match = re.match(r'([A-Z]+)(?=[a-z0-9_]|$)', name)
            if match:
                seq = match.group(1)
                if len(seq) > 1 and len(seq) < len(name): # MyURLIsNot -> myURLIsNot (preserve initialism)
                     return leading_underscores + seq[:-1].lower() + seq[-1] + name[len(seq):]
                elif len(seq) == len(name): # URL -> url, ABCD -> abcd
                     return leading_underscores + name.lower()
                else: # Myvar -> myvar (first letter cap, rest lower)
                     return leading_underscores + name[0].lower() + name[1:]
            else: # Should not happen if name[0].isupper(), but as fallback:
                return leading_underscores + name[0].lower() + name[1:]
        # if all lower (e.g. test), it's already camelCase by this definition
        return leading_underscores + name

    # Has underscores
    parts = name.split('_')
    parts = [part for part in parts if part]
    if not parts: return leading_underscores

    # Convert parts like "URL" in "my_URL_var" to "URL" not "Url" for capitalization phase
    # then the first part is lowercased. my_URL_var -> my + URL.capitalize() + Var.capitalize() -> myURLVar
    processed_parts = []
    for i, part in enumerate(parts):
        if part.isupper() and len(part) > 1: # like URL
            processed_parts.append(part)
        else:
            processed_parts.append(part.capitalize())

    return leading_underscores + processed_parts[0].lower() + ''.join(processed_parts[1:])


def to_pascal_case(name):
    if not name or name == "_" or (name.startswith("__") and name.endswith("__")):
        return name

    leading_underscores = ""
    idx = 0
    while idx < len(name) and name[idx] == '_':
        leading_underscores += "_"
        idx += 1
    name = name.lstrip('_')
    if not name: return leading_underscores

    if '_' not in name:
        if name[0].isupper():
            # If it's like "MyVAR" or "MYVar", we want "MyVar"
            # If "URL", keep "URL". If "MyURL", keep "MyURL".
            # This means if it starts with a cap, and isn't all caps, it's likely fine.
            # If it is all caps (URL), also fine.
            # The main conversion is for camelCase or all_lower.
             match = re.match(r'([A-Z]+)(?=[a-z0-9_]|$)', name)
             if match:
                seq = match.group(1)
                if len(seq) > 1 and len(seq) < len(name) : # MyURLIsNot -> MyURLIsNot (preserve initialism)
                    # Check if next char is upper, e.g. MyXYZVar -> MyXyzVar
                    # This is too complex for now. Let's assume MyURLIsNot is fine.
                    return leading_underscores + name
                # if len(seq) == len(name) -> URL, ABCD - keep as is for Pascal
                # if len(seq) == 1 -> Myvar - keep as is for Pascal
             return leading_underscores + name # if starts with Cap, assume ok or initialism.

        # if all lower (e.g. test -> Test) or camelCase (myVar -> MyVar)
        # For "url" -> "Url", "myUrl" -> "MyUrl"
        # For "myVAR" -> "MyVAR" (capitalize first letter)
        return leading_underscores + name[0].capitalize() + name[1:]

    # Has underscores
    parts = name.split('_')
    parts = [part for part in parts if part]
    if not parts: return leading_underscores

    # Convert parts like "URL" in "my_URL_var" to "URL" not "Url"
    # my_url_var -> MyUrlVar
    # my_URL_var -> MyURLVar
    pascal_parts = []
    for part in parts:
        if part.isupper() and len(part) > 1: # like URL
            pascal_parts.append(part)
        else:
            pascal_parts.append(part.capitalize())

    return leading_underscores + ''.join(pascal_parts)


class Scope:
    def __init__(self, parent=None, is_class_scope=False, class_name=None):
        self.parent = parent
        self.names_map = {}  # old_name -> new_name
        self.definitions = set() # names defined in this exact scope
        self.is_class_scope = is_class_scope
        self.class_name = class_name

    def define(self, old_name, new_name):
        self.names_map[old_name] = new_name
        self.definitions.add(old_name)

    def find_new_name(self, old_name):
        if old_name in self.names_map:
            return self.names_map[old_name]
        if self.parent:
            return self.parent.find_new_name(old_name)
        return None

    def is_defined_locally(self, old_name):
        return old_name in self.definitions

class NameTransformer(ast.NodeTransformer):
    def __init__(self):
        self.global_scope = Scope()
        self.current_scope = self.global_scope
        self.imported_names = set()
        self.class_attributes_map = {}
        self.renamed_globals = {} # old_global_name -> new_global_name (functions and classes)

    def visit_Import(self, node):
        for alias in node.names:
            name_to_add = alias.asname if alias.asname else alias.name.split('.')[0]
            self.imported_names.add(name_to_add)
        return self.generic_visit(node)

    def visit_ImportFrom(self, node):
        for alias in node.names:
            name_to_add = alias.asname if alias.asname else alias.name
            self.imported_names.add(name_to_add)
        return self.generic_visit(node)

    def enter_scope(self, is_class_scope=False, class_name=None):
        self.current_scope = Scope(parent=self.current_scope, is_class_scope=is_class_scope, class_name=class_name)

    def exit_scope(self):
        if self.current_scope.parent:
            self.current_scope = self.current_scope.parent

    def visit_ClassDef(self, node):
        original_name = node.name
        new_name = original_name

        if original_name not in self.imported_names and not (original_name.startswith('_') and not (original_name.startswith('__') and original_name.endswith('__'))): # Allow __ClassName__ but not _ClassName
            new_name = to_pascal_case(original_name)
            if new_name != original_name:
                self.current_scope.define(original_name, new_name)
                node.name = new_name
                if self.current_scope == self.global_scope:
                     self.renamed_globals[original_name] = new_name

        self.enter_scope(is_class_scope=True, class_name=new_name)
        if new_name not in self.class_attributes_map: # new_name is the potentially changed name
            self.class_attributes_map[new_name] = {}

        new_bases = [self.visit(base) for base in node.bases]
        node.bases = new_bases

        new_keywords = [ast.keyword(arg=kw.arg, value=self.visit(kw.value)) for kw in node.keywords]
        node.keywords = new_keywords

        new_decorator_list = [self.visit(d) for d in node.decorator_list]
        node.decorator_list = new_decorator_list

        new_body = []
        for stmt in node.body:
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)):
                new_body.append(self.visit_function_def_or_method(stmt, is_method=True, class_name=new_name))
            elif isinstance(stmt, ast.Assign):
                self.visit(stmt.value)
                for target in stmt.targets:
                    if isinstance(target, ast.Name):
                        class_var_original_name = target.id
                        if not class_var_original_name.startswith('_'):
                            class_var_new_name = to_camel_case(class_var_original_name)
                            if class_var_new_name != class_var_original_name:
                                self.current_scope.define(class_var_original_name, class_var_new_name) # Define in class scope
                                self.class_attributes_map[new_name][class_var_original_name] = class_var_new_name
                                target.id = class_var_new_name
                new_body.append(stmt)
            else:
                new_body.append(self.visit(stmt))
        node.body = new_body

        self.exit_scope()
        return node

    def visit_function_def_or_method(self, node, is_method=False, class_name=None):
        original_name = node.name
        new_name = original_name

        if not (original_name.startswith('__') and original_name.endswith('__')) and \
           not (is_method and original_name.startswith('_')) and \
           not (not is_method and original_name.startswith('_')) and \
           original_name not in self.imported_names:

            new_name = to_snake_case(original_name)
            if new_name != original_name:
                self.current_scope.define(original_name, new_name)
                node.name = new_name
                if self.current_scope == self.global_scope:
                    self.renamed_globals[original_name] = new_name
                elif is_method and class_name:
                    if class_name in self.class_attributes_map:
                         self.class_attributes_map[class_name][original_name] = new_name

        self.enter_scope(class_name=class_name if is_method else None)
        new_decorator_list = [self.visit(d) for d in node.decorator_list]
        node.decorator_list = new_decorator_list

        if node.args: self.visit(node.args)

        new_body = [self.visit(stmt) for stmt in node.body]
        node.body = new_body

        if node.returns: node.returns = self.visit(node.returns)
        self.exit_scope()
        return node

    def visit_FunctionDef(self, node):
        return self.visit_function_def_or_method(node, is_method=False)

    def visit_AsyncFunctionDef(self, node):
        return self.visit_function_def_or_method(node, is_method=False) # Assuming async functions can also be methods

    def visit_arguments(self, node):
        for arg_type_list in [node.posonlyargs, node.args, node.kwonlyargs]:
            for arg in arg_type_list:
                original_name = arg.arg
                if original_name.lower() not in ['self', 'cls'] and not original_name.startswith('_') and original_name not in self.imported_names:
                    new_name = to_camel_case(original_name)
                    if new_name != original_name:
                        self.current_scope.define(original_name, new_name)
                        arg.arg = new_name
                if arg.annotation: arg.annotation = self.visit(arg.annotation)

        for vararg_node_name in ['vararg', 'kwarg']:
            vararg_node = getattr(node, vararg_node_name)
            if vararg_node:
                original_name = vararg_node.arg
                if not original_name.startswith('_') and original_name not in self.imported_names:
                    new_name = to_camel_case(original_name)
                    if new_name != original_name:
                        self.current_scope.define(original_name, new_name)
                        vararg_node.arg = new_name
                if vararg_node.annotation: vararg_node.annotation = self.visit(vararg_node.annotation)

        for default_expr in node.defaults: self.visit(default_expr)
        for default_expr in node.kw_defaults:
            if default_expr: self.visit(default_expr)
        return node

    def visit_Name(self, node):
        original_id = node.id
        if original_id in self.imported_names: return node

        new_id_from_scope = self.current_scope.find_new_name(original_id)

        if new_id_from_scope:
            node.id = new_id_from_scope
        elif isinstance(node.ctx, ast.Store):
            if not original_id.startswith('_') and \
               not (original_id[0].isupper() and any(c.islower() for c in original_id[1:])) and \
               not original_id.isupper(): # Avoid renaming constants like ALL_CAPS if they are assigned.
                new_id_for_var = to_camel_case(original_id)
                if new_id_for_var != original_id:
                    self.current_scope.define(original_id, new_id_for_var)
                    node.id = new_id_for_var
        elif isinstance(node.ctx, ast.Load) and original_id in self.renamed_globals and not self.current_scope.find_new_name(original_id): # Check it's not shadowed
             node.id = self.renamed_globals[original_id]
        return node

    def visit_Attribute(self, node):
        node.value = self.visit(node.value)
        attr_name = node.attr

        # This part is heuristic, relying on `self` or direct ClassName.
        class_name_of_obj = None
        if isinstance(node.value, ast.Name):
            obj_name = node.value.id # This is the (potentially refactored) name of the object/class
            if self.current_scope.is_class_scope and obj_name == 'self':
                class_name_of_obj = self.current_scope.class_name
            elif obj_name in self.class_attributes_map: # Direct Class.attr access, obj_name is a known class name
                 class_name_of_obj = obj_name
            # TODO: Infer class_name_of_obj if obj_name is an instance variable. This needs type tracking.

        if class_name_of_obj and class_name_of_obj in self.class_attributes_map:
            if attr_name in self.class_attributes_map[class_name_of_obj]:
                node.attr = self.class_attributes_map[class_name_of_obj][attr_name]
            # If attr_name itself was a class variable that got camelCased
            # (less common for attributes to be non-snake/camel already, but possible)
            # This is mostly for methods that became snake_case.
        return node

    def visit_Assign(self, node):
        node.value = self.visit(node.value)
        new_targets = [self.visit(target) for target in node.targets]
        node.targets = new_targets
        return node

    def visit_Call(self, node):
        node.func = self.visit(node.func)
        new_args = [self.visit(arg) for arg in node.args]
        node.args = new_args
        new_keywords = [ast.keyword(arg=kw.arg, value=self.visit(kw.value)) for kw in node.keywords]
        node.keywords = new_keywords
        return node

def refactor_file_content(content, filepath="<string>"):
    try:
        tree = ast.parse(content, filename=filepath)
        transformer = NameTransformer()
        new_tree = transformer.visit(tree)
        ast.fix_missing_locations(new_tree)
        return ast.unparse(new_tree)
    except SyntaxError as e:
        return f"# Original file had SyntaxError: {e}\n{content}"
    except Exception as e:
        # import traceback
        # tb_str = traceback.format_exc() # Keep this for debugging if the agent fails
        return f"# Error during refactoring of {filepath}: {e}\n# AST processing failed.\n{content}"

def get_python_files(directory):
    py_files = []
    for root, dirs, files in os.walk(directory):
        # Efficiently remove __pycache__ from traversal
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for file in files:
            if file.endswith(".py"):
                py_files.append(os.path.join(root, file))
    return py_files
