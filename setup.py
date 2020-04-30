from setuptools import setup, find_packages
import versioneer

requirements = [
    "pyyaml",
    "pyspark",
    "glow.py",
]

setup(
    name='firefly',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Short description",
    license="MIT",
    author="Florian R. Hölzlwimmer",
    author_email='git.ich@frhoelzlwimmer.de',
    url='https://github.com/hoeze/firefly',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'firefly=firefly.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='firefly',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ]
)
