from setuptools import setup
import os

if __name__ == "__main__":
    os.system("pip install numpy==1.19.3 --only-binary=:all:")
    setup(include_package_data=True)
