import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DARt-Yuting-xxxxxx", # Replace with your own username
    version="0.0.1",
    author="Xian Yang, Shuo Wang, Yuting Xing, Yida Xu, Yike Guo",
    author_email="author@example.com",
    description="A tool to estimate Instantaneous Reproduction Number(Rt) for the pandemic",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)