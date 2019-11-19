# Differential_Expression-Aging_Brain_Tissue
Time is the grim-reaper to us all. As we age, our cells accumulate mutations, oxidative damage, and (other examples - general entropy) which  affects our cells to function properly, cascading to tissue and organ damage. Using RNA-seq data from the open source Genome-Tissue Expression project at the Broad Institute, we are interested in measuring differential expression of genes (will find subset for this project) in brain tissues amongst aging individuals.
(This is just a start, will look over articles for more cohesive abstract)
## Getting Started

- A de-identified, open access version of the sample annotations available in dbGaP. We use this file to extract the brain tissue types of individuals. (11M)
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
```
    - This file was further simplified by extracting the desired tissue types (all brain tissue and nerve tissue). The resulting file (ID_Brain_Nerve_Tissue.txt) is used in the main python code. 
    ```
    $ awk '/Brain/ || /Nerve/' GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | awk -F'\t' '{print $1, $7}' OFS='\t'  > ID_Brain_Nerve_Tissue.txt
    ```

- A de-identified, open access version of the subject phenotypes available in dbGaP. We use this file to extract the ages of individuals. (20K)
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
```

- Gene TPMs. We use this file to extract the gene transcript per million (TPM) for an individual.
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
```

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
