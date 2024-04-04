from setuptools import setup, find_packages

setup(name='SiteSeeker',
    version='1.0',
    description='Machine-learning structure-based predictron of protein ligand binding sites',
    author='Arnau Noguera Segura, Aina Vaquer Picó',
    author_email='arnau.noguera01@estudiant.upf.edu, aina.vaquer01@estudiant.upf.edu',
    url='https://github.com/arnaunoguera/SiteSeeker',
    install_requires=['biopython>=1.83','numpy>=1.26.4','pandas>=2.2.1','joblib>=1.3.2','scikit-learn>=1.4.1.post1'],
    extras_require={'dssp': ['dssp']},
    packages=['SiteSeeker', 'Training'],
    package_data={'SiteSeeker': ['model.joblib', 'scaler.pkl'],
                  'Training': ['hyperparameter_tuning.tsv', 'ligands.txt', 'training_file.tsv']},
    entry_points={
        'console_scripts': [
            'SiteSeeker=SiteSeeker.PredictBindingSites:main'
        ]
    }
)