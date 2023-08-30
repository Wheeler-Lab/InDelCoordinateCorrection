import pathlib

# Class that takes the config and processes the tree of experiments
# for later use in the snakemake rules
class ExperimentTree:
    def __init__(self, config):
        # Results path
        self.results = pathlib.Path('results')
        # Samples directory
        self.samples_directory = pathlib.Path(config['samples_directory'])
        # Experiments specification in the configuration
        self.experiments: dict = config['experiments']
        # Experiments to build (everything that has a parent)
        self.to_build: list[str] = [k for k, v in self.experiments.items() if 'parent' in v]

        # Process the dependency tree
        self._digest_dependencies()
        # Create a lookup table to find which fasta file goes with which sequencing reads
        self._digest_fasta_lookup()
        # Create the directory structure
        self._make_directory_structure()

    def _digest_dependencies(self):
        self.dependencies: dict = {
            experiment: {
                "samples": self.experiments[experiment]['samples']
            }
            for experiment in self.to_build
        }

        for experiment, entry in self.dependencies.items():
            parent_id = self.experiments[experiment]['parent']
            parent = self.experiments[parent_id]
            if 'parent' not in parent:
                entry['input_gff'] = parent['genome_annotations']
                entry['input_fasta'] = parent['genome_sequence']
            else:
                entry['input_gff'] = self.results / parent_id / 'genome.gff'
                entry['input_fasta'] = self.results / parent_id / 'genome.fasta'

    def _digest_fasta_lookup(self):
        self.fasta_lookup: dict[str, str] = {
            str(self.results / experiment / sample): self.dependencies[experiment]['input_fasta']
            for experiment in self.to_build for sample in self.dependencies[experiment]['samples']
        }

    def _make_directory_structure(self):
        for experiment in self.to_build:
            (self.results / experiment).mkdir(parents=True, exist_ok=True)

    @property
    def all_inputs(self):
        return [self.results / experiment / "genome.gff" for experiment in self.to_build]

