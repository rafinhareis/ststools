from setuptools import setup

with open("README.md", "r") as arq:
    readme = arq.read()

setup(name='ststools',
    version='1.0.1',
    license='MIT License',
    author='Rafael Reis Barreto',
    long_description=readme,
    long_description_content_type="text/markdown",
    author_email='rafinhareis17@gmail.com',
    keywords='nanosurf,omicron,sts',
    description=u'Visualizacao de sts nanosurf e omicron',
    packages=['ststools'],
    install_requires=['pandas','numpy','scipy','matplotlib','ipywidgets','Pillow','scikit-learn'],)