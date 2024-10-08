import ast
import os


PACKAGE = 'apero'
CHAINS = set()

def find_imports(file_path):
    """Parse a Python file and return all imported modules."""
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    imports = set()

    for line in lines:
        if line.startswith('import '):
            module = line.split('import ')[1].split()[0]

            if PACKAGE is not None and PACKAGE in module:
                imports.add(module)
        elif line.startswith('from '):
            module1 = line.split('from ')[1].split()[0]
            module2 = line.split('import ')[1].split()[0]
            module = f"{module1}.{module2}"
            if PACKAGE is not None and PACKAGE in module:
                imports.add(module)

    return imports


def trace_dependencies(file_path, base_path, tree=None, circular=None, chain=None):
    """Trace dependencies hierarchically and handle circular imports."""

    if tree is None:
        tree = {}
    if circular is None:
        circular = set()
    if chain is None:
        chain = []

    file_name = os.path.relpath(file_path, base_path)

    # Check for circular dependency
    if file_name in chain:
        # Mark circular dependencies in the tree
        tree[file_name] = "Circular"
        circular.add(file_name)
        return tree, circular

    # Add to current chain
    chain.append(file_name)

    if file_name not in tree:
        tree[file_name] = {}

    imports = find_imports(file_path)

    for module in imports:
        module_path = os.path.join(base_path, f"{module.replace('.', os.sep)}.py")
        module_path2 = os.path.dirname(module_path) + '.py'
        module_path3 = os.path.dirname(module_path2) + '.py'
        # Check if the module is a python file
        if os.path.exists(module_path):
            tree[file_name][module], circ = trace_dependencies(module_path, base_path, {}, circular, chain.copy())
            circular = circular | circ
        # Check if level up is the python file (and isn't the same as file_name)
        elif os.path.exists(module_path2) and module_path2 != file_path:
            tree[file_name][module], circ = trace_dependencies(module_path2, base_path, {}, circular, chain.copy())
            circular = circular | circ
         # Check if level up
        elif os.path.exists(module_path3) and module_path3 != file_path:
            tree[file_name][module], circ = trace_dependencies(module_path3, base_path, {}, circular, chain.copy())
            circular = circular | circ
        else:
            # If the file doesn't exist, just add the module name
            tree[file_name][module] = {}

    return tree, circular


def print_dependency_tree(all_lines, tree, indent=0, prefix="", display_connections=True):
    """Helper function to print the hierarchical dependency tree with guides."""

    if display_connections:
        c1 = "|___ "
        c2 = "|--- "
        c3 = "|   "
    else:
        c1 = "     "
        c2 = "     "
        c3 = "    "

    extension = ''
    for i, (key, value) in enumerate(tree.items()):
        connector = c1 if i == len(tree) - 1 else c2
        all_lines += [f"{prefix}{connector}{key}"]

        if isinstance(value, dict):
            # Continue with deeper levels
            extension = "    " if i == len(tree) - 1 else c3
            print_dependency_tree(all_lines, value, indent + 1,
                                  prefix + extension,
                                  display_connections=display_connections)
        else:
            all_lines += [f"{prefix}{extension}{value}"]
    return all_lines


def module_dependencies(code_file, base_path=None):
    """Main function to find all dependencies of a module."""
    if base_path is None:
        base_path = os.path.dirname(code_file)
    tree, circular = trace_dependencies(code_file, base_path)

    print("Dependency Tree:")
    all_lines1 = print_dependency_tree([], tree)
    # write all lines to text file
    with open('dependency_tree.txt', 'w') as f:
        for line in all_lines1:
            f.write(line + '\n')
    all_lines2 = print_dependency_tree([], tree, display_connections=False)
    # write all lines to text file
    with open('dependency_tree_clean.txt', 'w') as f:
        for line in all_lines2:
            f.write(line + '\n')


    if circular:
        print("\nCircular dependencies detected in the following paths:")
        for cycle in circular:
            print(" -> ".join(cycle))
    else:
        print("\nNo circular dependencies detected.")



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":



    module_dependencies('/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/'
                        'apero/tools/recipes/bin/apero_validate.py',
                        '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/')


