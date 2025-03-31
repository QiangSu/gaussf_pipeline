from setuptools import setup, find_packages

setup(
    name="gaussf_pipeline",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=["pandas", "scipy", "numpy", "tqdm"],
)
