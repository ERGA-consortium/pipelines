# ERGA AUGUSTUS resources

A collection of AUGUSTUS model parameters trained by members of the ERGA Annotation Community. 
If you wish to contribute a parameter file, please fork the pipelines repository and create a pull request.
The structure of this folder is as follows:

```
|-AUGUSTUS
| |-species_name
| | |-README.md
| | |-AUGUSTUS.params
```

Please ensure to document how you created your training file, in particular:
* Which assembly was used, which GCA accession number
* Which evidences were used, e.g. RNA-seq, Iso-seq, protein sequences with accession numbers
* Which tools were used to train the model, with version numbers
* Links to any publications associated with the model
* A list of authors and contact details should quesions arise in the future
