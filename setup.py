from setuptools import setup

setup(name='penta',
      version='0.1',
      description='a toy HF code',
      url='http://github.com/ZhaoYilin/penta',
      author='Yilin Zhao',
      author_email='zhaoyilin1991@gmail.com',
      license='MIT',
      packages=['penta'],
      test_suite='nose.collector',
      install_requires=[
        'numpy',
        'nose',
      ],
      zip_safe=False)
