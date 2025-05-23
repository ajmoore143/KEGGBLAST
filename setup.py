#!/usr/bin/env python
# coding: utf-8

# In[1]:


from setuptools import setup, find_packages

setup(
    name="keggblast",
    version="0.1",
    description="KEGG Orthology Gene Extraction + BLAST Toolkit",
    author="Alex Moore",
    packages=find_packages(),
    install_requires=[
        "requests",
        "pandas",
        "gget",
        "rapidfuzz"
    ],
    python_requires=">=3.7"
)


# In[ ]:




