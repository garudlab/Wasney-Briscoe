from setuptools import setup, find_packages

# List of requirements
requirements = []  # This could be retrieved from requirements.txt
# Package (minimal) configuration
setup(
    name="PopGenCalculator",
    version="1.0.0",
    description="devo",
    packages=find_packages(),  # __init__.py folders search
    install_requires=requirements,
)
