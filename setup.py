from setuptools import setup,find_packages
setup(
    name="Utility",
    version="1.2.3",
    author="Shinsuke Sakai",
    url="https://github/ShinsukeSakai0321/Utility.git",
    packages=find_packages("src"),
    package_dir={"":"src"},
    package_data={'':['*.csv']},
    include_package_data=True,
)