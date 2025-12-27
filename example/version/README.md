# Version Example

This folder contains a simple example that demonstrates how to query the version information of the `t8code` library that your application is linked against.

## Files

- `t8_version.cxx`: This program calls the `t8code` version functions and prints the library's version number, git commit hash, and other relevant information to the console.

## Main Purpose

This example showcases a simple but important piece of functionality. Being able to programmatically access the library version is crucial for:
- **Reproducibility:** Logging the exact version of `t8code` used in a simulation is essential for ensuring that results can be reproduced later.
- **Debugging:** When reporting bugs, knowing the precise version of the library is vital for developers to understand the context of the issue.
- **Feature Checking:** An application could potentially check the `t8code` version to ensure that a required feature is available.

## How to Run

After building the examples, you can simply run the executable:

```bash
./t8_version
```

It will print the `t8code` version details to standard output.
