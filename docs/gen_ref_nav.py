"""Generate the code reference pages and navigation.
Assumptions:
- Code that we want to autoref is accessible/importable by pytkdoc. In practice, this means that the module is installed, probably via a `pip editale`
- The documentation pages will be put under `ref_dir`
- The source code file structure is:

```
|-- docs
|   |-- index.md
|   |-- reference/
|-- src
    |-- code.py
    |-- package
        |-- morecode.py
```

"""
### ONLY USER SETTINGS:
srcdir = 'src'
ref_dir = 'reference'
ignore_files = ['setup.py','__']
"""
skeleton borrowed from https://github.com/mkdocstrings/mkdocstrings/blob/master/docs/gen_ref_nav.py
"""

### REUSABLE CODE:
from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

for path in sorted(Path(srcdir).glob("**/*.py")): #search over all subdirs for all .py's
    ignore = False
    for ignore_this in ignore_files: 
        if ignore_this in str(path):
            ignore = True
            break
    if ignore: continue

    module_path = path.relative_to(srcdir).with_suffix("") #need to go *up* one, assuming that the setup.py is in...
    doc_path = path.relative_to(srcdir).with_suffix(".md")
    full_doc_path = Path(ref_dir, doc_path)

    #print(module_path)
    #print(doc_path)
    #print(full_doc_path)

    parts = list(module_path.parts)
    parts[-1] = f"{parts[-1]}.py"
    nav[parts] = doc_path

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        identifier = ".".join(module_path.parts)
        print("::: " + identifier, file=fd)

    mkdocs_gen_files.set_edit_path(full_doc_path, path)

with mkdocs_gen_files.open("{}/SUMMARY.md".format(ref_dir), "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())

''' Not used:
nav["mkdocs_autorefs", "references.py"] = "autorefs/references.md"
nav["mkdocs_autorefs", "plugin.py"] = "autorefs/plugin.md"
'''
