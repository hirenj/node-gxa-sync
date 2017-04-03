"use strict";
// ES6

const request = require('request');

const experiments_url = 'http://www.ebi.ac.uk/gxa/json/experiments';

const ftp = require('ftp');
const parse_url = require('url').parse;
const csv = require('csv-parse');
const xmlNodes = require('xml-nodes');
const xmlObjects = require('xml-objects');

const baseline_whitelist = [ 'E-MTAB-5214', 'E-MTAB-2836'];

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
    return [].concat((exp.species || [])).filter((tax) => species_set.has(tax)).length > 0;
  });
};

let filterBaseline = (experiments) => {
  return experiments.filter((exp) => {
    return exp.experimentType === 'RNASEQ_MRNA_BASELINE';
  }).filter( exp => {
    return (baseline_whitelist.length == 0) || baseline_whitelist.indexOf(exp.experimentAccession) >= 0;
  });
};

let getExperimentInfo = (exp) => {
  return {
    accession: exp.experimentAccession,
    date: exp.lastUpdate,
    sources: exp.experimentalFactors,
    taxonomy: [].concat(exp.species || []).join(','),
    description: exp.experimentDescription
  };
};

let getBaselineForSpecies = (species) => {
  return getExperiments()
  .then(filterSpecies.bind(null,species))
  .then(filterBaseline)
  .then((exps) => exps.map(getExperimentInfo));
};

let make_ontology_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}.condensed-sdrf.tsv`;
};

let make_configuration_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}-configuration.xml`;
};

let make_data_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}.tsv`;
};

const read_configuration = function(stream) {
  let result = {};
  let data_stream = stream.pipe(xmlNodes('assay_group'))
  .pipe(xmlObjects({explicitRoot: false, explicitArray: false, mergeAttrs: true}))
  .on('data', dat => {
    result[dat.id] = { label: dat.label, sample: dat.assay[0]['_'] };
  });
  return new Promise( (resolve,reject) => {
    data_stream.on('end',resolve);
    stream.on('error',reject);
    data_stream.on('error',reject);
  }).then( () => result );
};


const read_ontology = function(stream) {
  return new Promise( (resolve,reject) => {
    let parser = csv({delimiter: "\t", relax_column_count: true}, function(err, data){
      if (err) {
        reject(err);
        return;
      }
      let ids = data.filter( row => row[3] == 'characteristic' && row[4] == 'organism part')
      .map( row => { return { ontology: row[6], sample_id: row[2]  } })
      .filter( purl => purl && purl.ontology )
      .map( purl => {
        purl.ontology = purl.ontology.split('/').reverse().shift().replace('_',':');
        return purl;
      });
      resolve(ids);
    });
    stream.pipe(parser);
  });
};

const read_data = function(stream) {
  return new Promise( (resolve,reject) => {
    resolve(stream.pipe(csv({delimiter: "\t", relax_column_count: true, columns: true})));
  });
};

const connect_ftp = function(url) {
  let req = parse_url(url);
  let ftp_site  = new ftp();
  let result = new Promise(function(resolve,reject) {
    ftp_site.on('ready',function() {
      resolve(ftp_site);
    });
    ftp_site.on('error',reject);
  });
  req.connTimeout = 1000;
  ftp_site.connect(req);
  return result.catch(function(err) {
    console.log(err);
    return null;
  });
};

const download_file = function(url,ftp) {
  let path = parse_url(url).path;
  return new Promise(function(resolve,reject) {
    ftp.get(path, function(err, stream) {
      if (err) {
        reject(err);
      } else {
        resolve(stream);
      }
    });
  });
};

let downloadExpData = function(experiment_id) {
  let ontology_url = make_ontology_url(experiment_id);
  console.log(ontology_url);
  let ontology_ids = connect_ftp(ontology_url)
  .then( download_file.bind(null,ontology_url) )
  .then( read_ontology )
  // .then( ids => console.log(ids) )
  .catch( console.log.bind(console) );

  // Download files to a local folder if they are newer than a given time..

  // ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-3358/E-MTAB-3358.tsv
  let data_url = make_data_url(experiment_id);
  console.log(data_url);
  let data = connect_ftp(data_url)
  .then( download_file.bind(null,data_url) )
  .then( read_data )
  .then( data_stream => data_stream.on('data', dat => console.log(dat) ));

/*
{ 'Gene ID': 'ENSG00000000938',
  'Gene Name': 'FGR',
  g6: '5,7,9,12,15',
  g10: '2,2,2,3,4',
  g12: '2,3,5,6,7',
  g7: '1,2,3,3,3',
  g18: '3,4,5,14,22',
  g3: '1,1,2,2,3',
  g26: '15,19,21,23,29',
  g19: '0.9,0.9,1,1,1',
  g30: '9,9,9,9,9',
  g27: '0.9,2,2,2,2',
  g20: '0.7,1,2,2,2',
  g4: '1,1,1,1,1',
  g1: '14,15,18,21,21',
  g8: '0.9,2,2,3,4',
  g28: '0.4,0.8,1,2,4',
  g25: '0.8,0.8,0.8,0.9,1',
  g13: '2,4,5,22,33',
  g2: '18,34,49,56,62',
  g23: '0.7,0.7,0.8,1,2',
  g15: '1,2,2,3,4',
  g17: '0.7,0.9,1,1,1',
  g31: '45,45,50,55,55',
  g24: '0.9,1,2,2,2',
  g21: '1,2,2,2,3',
  g9: '3,3,3,4,5',
  g16: '52,55,101,145,145',
  g32: '3,3,3,4,7',
  g29: '2,2,2,2,2',
  g14: '1,3,3,4,6',
  g11: '2,3,4,4,4',
  g22: '0.2,0.2,0.3,0.5,0.7',
  g5: '11,11,12,13,14' }
*/

  // Entries have min,percentile,median,percentile,max for each row apparently?

  let configuration_url = make_configuration_url(experiment_id);
  console.log(configuration_url);
  let configuration = connect_ftp(configuration_url)
  .then( download_file.bind(null,configuration_url) )
  .then( read_configuration )
  .then( groups => console.log(groups) )

  // Parse the files, extracting the info we want and translating the ensembl to entrez ids

  // Write this file out to a summarised data folder (in { data: {}, metadata: {} } format)

  // Let Gator daemon figure out the hosting duties
};

let foo = getBaselineForSpecies(['Homo sapiens']).then((exps) => {
  console.log(exps);
  downloadExpData(exps[0].accession);
  return exps;
}).catch((err) => console.error(err));