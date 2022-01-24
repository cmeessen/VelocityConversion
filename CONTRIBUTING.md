# Contributing

To contribute to `VelocityConversion`, feel free to open an
[issue](https://github.com/cmeessen/VelocityConversion/issues) anytime. Of
course, you are also welcome to open pull requests from your personal forks.

## Getting started

First, create your own personal fork of the repository using the fork button in
the upper right corner of this page. Now, clone the fork to your personal hard
drive:

```bash
# Note, that the URL will differ for your repository
git clone git@github.com:cmeessen/VelocityConversion.git
```

or via HTTPS:

```bash
git clone https://github.com/cmeessen/VelocityConversion.git
```

and `cd` into the corresponding folder on your disk

```bash
cd VelocityConversion
```

### Create a virtual environment

Use pipenv to create a new virtual environment

```bash
pipenv shell
```

and install the dev-dependencies:

```bash
pipenv shell
pipenv install --dev
```

## Code style

The python code in this repository is checked with `pycodestyle`. Please The rules
for `pycodestyle` are defined in `setup.cfg`. Before creating a pull request,
please check your code by running:

```bash
pycodestyle
```

If the command does not generate any output, the style is correct.

## Building documentation

To build the documentation, go to the `docs/sphinx` folder,

```bash
cd docs/sphinx
```

and run

```bash
make html
```
