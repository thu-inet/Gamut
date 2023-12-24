'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2023-12-24 11:51:31
Description: 

Copyright (c) 2023 by THU-RSAG, All Rights Reserved. 
'''
import setuptools

with open("README.md", "r") as fileopen:
    long_description = fileopen.read()

setuptools.setup(
    name="Gamut",
    version="0.1.1",
    author="THU-RSAG",
    author_email="z-wj21@mails.tsinghua.edu.cn",
    description="Gamut is a Python package for the analysis of 1-d gamma spectrum.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/thu-inet/Gamut",
    packages=setuptools.find_packages(),
    install_requires=['numpy>=1.25.0', 'matplotlib>=3.7.0'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
)
