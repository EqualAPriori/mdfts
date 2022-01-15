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
src_dir = 'src'
ref_dir = 'reference'
ignore_files = ['setup.py','__']
ignore_dirs = ['.egg-info','__']
"""
skeleton borrowed from https://github.com/mkdocstrings/mkdocstrings/blob/master/docs/gen_ref_nav.py
"""

### REUSABLE CODE:
from pathlib import Path

import os
import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

# FOR CUSTOM ToC
paths = []
levels = []
relpaths = []
links = []
indentation = 4

# ITERATE
for path in sorted(Path(src_dir).glob("**/*")): 
    #if non-ignored dir, save
    if path.is_dir():
        if bool( [ele for ele in ignore_dirs if (ele in path.name)] ):
            continue
        else: 
            level = str(path.relative_to(src_dir)).count(os.sep)
            levels.append(level)
            relpaths.append(path.name)
            links.append('')
            continue
    #check if is a file that we want to ignore 
    if not path.name.endswith('py'):
        continue
    if bool( [ele for ele in ignore_files if (ele in path.name)] ): #check if elements that we want to ignore occur
        continue

    paths.append(str(path))
    module_path = path.relative_to(src_dir).with_suffix("") 
    doc_path = path.relative_to(src_dir).with_suffix(".md")
    full_doc_path = Path(ref_dir, doc_path)

    level = str(path.relative_to(src_dir)).count(os.sep)
    levels.append(level)
    relpaths.append(path.name)
    links.append(doc_path)


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
    #print('\n'.join(list(nav.build_literate_nav())))

with mkdocs_gen_files.open("{}/ToC.md".format(ref_dir),'w') as toc_file:
    print("===")
    print('AUTOREF ToC: ')
    toc_file.write("# Code Tree\n")
    for lvl,relpath,link in zip( levels,relpaths,links ):
        if link == "":
            line = '{sep}- {name}/'.format(sep=lvl*' '*indentation,name=relpath)
        else:
            line = '{sep}- [{name}]({link})'.format(sep=lvl*' '*indentation,name=relpath,link=link)
        print(line)
        toc_file.write(line+"\n")
    print("===")
    #toc_file.writelines("\n".join(paths))

''' Not used:
nav["mkdocs_autorefs", "references.py"] = "autorefs/references.md"
nav["mkdocs_autorefs", "plugin.py"] = "autorefs/plugin.md"
'''
