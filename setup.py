from setuptools import find_packages, setup
import os

# Set the base directory location which is the location of the setup.py file
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

setup(
    name="phanta",
    version="0.1",
    description="Python package to assemble phage illumina reads",
    author="Joshua J Iszatt",
    author_email="joshiszatt@gmail.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
    ],
    python_requires=">=3.10",
    packages=find_packages(include=['phanta', 'phanta.*']),
    package_data={
        'phanta': ['config.json'],
    },
)

