"use strict";
// ES6

const request = require('request');

const experiments_url = 'http://www.ebi.ac.uk/gxa/json/experiments';

let getExperiments = () => {
	return new Promise(function(resolve,reject) {
		request(experiments_url,function(err,response,body) {
			if ( err ) {
				return resolve(err);
			}
			let parsed = JSON.parse(body);
			if ( ! parsed.aaData) {
				return resolve(new Error('Invalid data format'));
			}
			resolve(parsed.aaData);
		});
	});
};

let filterSpecies = (species,experiments) => {
	let species_set = new Set(species);
	return experiments.filter((exp) => {
		return (exp.species || []).filter((tax) => species_set.has(tax)).length > 0;
	});
};

let filterBaseline = (experiments) => {
	return experiments.filter((exp) => {
		return exp.experimentType === 'RNASEQ_MRNA_BASELINE';
	});
};

let getExperimentInfo = (exp) => {
	return {
		accession: exp.experimentAccession,
		date: exp.lastUpdate,
		sources: exp.experimentalFactors,
		taxonomy: exp.species.join(','),
		description: exp.experimentDescription
	};
};

let getBaselineForSpecies = (species) => {
	return getExperiments().then(filterSpecies.bind(null,species)).then(filterBaseline).then((exps) => exps.map(getExperimentInfo));
};

getBaselineForSpecies = function() {
	return Promise.resolve(require('./exps.parsed.json'));
};

downloadExpData = function() {
	// Possibly not needed if we just trust the text for the organism part..
	// ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-3358/E-MTAB-3358.condensed-sdrf.tsv

	// Download files to a local folder if they are newer than a given time..

	// ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-3358/E-MTAB-3358.tsv

	// Entries have min,percentile,median,percentile,max for each row apparently?

	// ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-3358/E-MTAB-3358-configuration.xml

	// Parse the files, extracting the info we want and translating the ensembl to entrez ids
	
	// Write this file out to a summarised data folder (in { data: {}, metadata: {} } format)

	// Let Gator daemon figure out the hosting duties
};

let foo = getBaselineForSpecies(['Homo sapiens']).then((exps) => {
	console.log(exps);
	return exps;
}).catch((err) => console.error(err));