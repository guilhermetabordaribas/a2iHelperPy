import setuptools

with open('README.md','r') as fh:
    description = fh.read()

setuptools.setup(
    name='GoTx',
    version='0.0.1',
    author='Guilherme Taborda Ribas',
    author_email='guilherme.ribas@gmail.com',
    packages=['gotx'],
    description=['A package to analyze microarray and RNA-seq and run machine learning algorithms.'],
    long_description=description,
    long_description_content_type='text/markdown',
    url='https://github.com/guilhermetabordaribas',
    license='MIT',
    python_requires='>=3.10',
    install_requires=[]
)
