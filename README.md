# qIT_Auto_Analyser

## Install

- To install directly from github use the following command:

```bash
pipx install --python python3.13 git+https://github.com/charlieb555/qIT_Auto_Analyser.git
```

This installs the application in the `~/.local/bin` and sets up the correct python virtual environment.  If `~/.local/bin` is not currently in your PATH environment variable, you need to add something like the following to the correct startup script for whatever shell you are using (e.g. `~/.bashrc` for bash, `~/.zshrc` for zsh):

```bash
export PATH="~/.local/bin:$PATH"
```

- You can check for updates with the following command:

```bash
pipx upgrade qIT_Auto_Analyser
```

- The following command can be used to later uninstall this application:

```bash
pipx uninstall qIT_Auto_Analyser
```
