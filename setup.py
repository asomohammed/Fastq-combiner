from setuptools import setup, find_packages

setup(
    name='fastq-combiner',
    version='1.0.0',
    description='High-speed, memory-efficient FASTQ file combiner for Cell Ranger and NGS pipelines',
    author='Your Name',
    packages=find_packages(),
    py_modules=['fastq_combiner'],
    install_requires=[
        'tqdm',
        'pyyaml',
    ],
    entry_points={
        'console_scripts': [
            'fastq-combiner=fastq_combiner:main',
        ],
    },
    python_requires='>=3.7',
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
) 