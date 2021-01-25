from setuptools import find_packages, setup


def read(path):
    with open(path, 'r') as f:
        return f.read()


setup(
    name='dynast',
    version='0.0.1',
    url='',
    author='Kyung Hoi (Joseph) Min',
    author_email='phoenixter96@gmail.com',
    description='',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    keywords='',
    python_requires='>=3.7',
    license='MIT',
    packages=find_packages(exclude=('tests', 'tests.*')),
    zip_safe=False,
    include_package_data=True,
    install_requires=read('requirements.txt').strip().split('\n'),
    entry_points={
        'console_scripts': ['dynast=dynast.main:main'],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
    ],
)
