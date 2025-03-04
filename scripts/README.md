# t8code scripts

This folder contains several scripts that are useful for t8code users and developers.

## Indentation

The purpose of the indentation scripts is to help t8code developers to indent their code according to the t8code [indentation guidelines](https://github.com/DLR-AMR/t8code/wiki/Coding-Guideline#indentation). Please read these guidelines before using these scripts.

#### t8indent

Apply this script to a `.c`, `.cxx`, `.h` or `.hxx` file to indent it according to the guidelines.
You can also provide a list of files and all will get indented.

This script uses the `clang-format` program with t8code specific settings.
Sometimes `t8indent` does produce undesired results. Therefore, after indenting use `git add -p` or similar to check all changes before committing. You can protect code lines from being changed by the script by enclosing them in `/* clang-format off */` and `/* clang-format on */` comments.

#### pre-commit

This script should be copied to your `.git/hooks` folder. `git` then automatically checks the indentation of committed files and prevents you from committing wrongly indented files. See [Git indentation workflow](https://github.com/DLR-AMR/t8code/wiki/Coding-Guideline#git-indentation-workflow).

#### check_if_file_indented.sh

Check for a single source file whether it is indented according to `t8indent`.

#### check_if_all_files_indented.sh

Check whether all t8code source files are properly indented.

#### indent_all_files.sh

This script indents all t8code source files at once. This script should only be used by the main developers of t8code. Handle with care.

#### find_all_source_files.sh

List all source files of t8code in the `src/` `example/` and `test/` subfolders.

## Others

#### su2_mesh_to_gmsh.py

Used to convert coarse meshes from an su2 format into msh format.

#### create_todo_list.sh

Parses through the code and collects comments containing the `TODO` keyword.
