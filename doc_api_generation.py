


import os
import ast

SRC_ROOT = "src/pyfastwell"
DOCS_ROOT = "docs/api"

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

all_items = []

for dirpath, _, filenames in os.walk(SRC_ROOT):
    rel_dir = os.path.relpath(dirpath, SRC_ROOT)
    docs_dir = os.path.join(DOCS_ROOT, rel_dir)
    ensure_dir(docs_dir)
    for filename in filenames:
        if filename.endswith(".py"):
            module_name = filename[:-3]
            file_path = os.path.join(dirpath, filename)
            with open(file_path, "r", encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=file_path)
            for node in tree.body:
                if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
                    name = node.name
                    if not name.startswith("_"):  # Skip private members and functions
                        md_path = os.path.join(docs_dir, f"{name}.md")
                        module_path = os.path.join(rel_dir, module_name).replace(os.sep, ".")
                        with open(md_path, "w", encoding="utf-8") as md:
                            md.write(f"# {name} API Reference\n\n")
                            md.write(f"::: pyfastwell.{module_path}.{name}\n")
                        rel_md_path = os.path.relpath(md_path, DOCS_ROOT)
                        all_items.append((name, rel_md_path))

# Create root API index
index_path = os.path.join(DOCS_ROOT, "index.md")
with open(index_path, "w", encoding="utf-8") as idx:
    idx.write("# API Reference Index\n\n")
    for name, rel_md_path in sorted(all_items):
        idx.write(f"- [{name}]({rel_md_path})\n")