# t8code scripts

## Indentation

The purpose of the indentation scripts is to help t8code developers to indent their code according to the t8code [indentation guidelines](https://github.com/holke/t8code/wiki/Coding-Guideline#indentation). Please read these guidelines before using these scripts.

#### t8indent

Apply this script to a `.c`, `.cxx`, `.h` or `.hxx` file to indent it according to the guidelines.
You can also provide a list of files and all will get indented.

This script uses the `GNU indent` program with t8code specific settings.
Sometimes `t8indent` does produce undesired results. Therefore, after indenting use `git add -p` or similar to check all changes before committing. You can protect code lines from being changed by the script by enclosing them in `/* *INDENT-ON* */` and `/* *INDENT-OFF* */` comments.
See also [Known issues with indent](https://github.com/holke/t8code/wiki/Known-issues-with-the-indent-script).

#### pre-commit

This script should be copied to your `.git/hooks` folder. `git` the automatically checks the indentation of committed files and prevents you from committing wrongly indented files. See [Git indentation workflow](https://github.com/holke/t8code/wiki/Coding-Guideline#git-indentation-workflow).

#### check_if_file_indented.scp

Check for a single source file whether it is indented according to `t8indent`.

#### check_if_all_files_indented.scp

Check whether all t8code source files are properly indented.

#### indent_all_files.scp

This script indents all t8code source files at once. This script should only be used by the core developers. Handle with care.

#### find_all_source_files.scp

List all source files of t8code in the `src/` `example/` and `test/` subfolders.

#### remove_double_const.scp

One of the [Known issues with indent](https://github.com/holke/t8code/wiki/Known-issues-with-the-indent-script) is that it duplicates the `const` keyword in `C++` files in some cases.

## Others

- su2_mesh_to_gmsh.py

- create_todo_list.sh
