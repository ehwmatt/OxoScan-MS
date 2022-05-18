import setuptools
import os

base_dir = os.path.dirname(os.path.realpath(__file__))

about = {}
with open(os.path.join(base_dir, "glycoproteomics", "__about__.py"), encoding="utf-8") as in_f:
    exec(in_f.read(), about)

setuptools.setup(
    name=about["__title__"],
    version=about["__version__"],
    author=about["__author__"],
    description=about["__summary__"],
    url=about["__url__"],
    package_dir={"": "."},
    packages=["glycoproteomics"],
    python_requires=">=3.9",
)
