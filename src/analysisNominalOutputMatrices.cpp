//FastQTL: Fast and efficient QTL mapper for molecular phenotypes
//Copyright (C) 2015 Olivier DELANEAU, Alfonso BUIL, Emmanouil DERMITZAKIS & Olivier DELANEAU
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "data.h"
#include "cnpy.h"

void data::runNominalOutputMatrices(string dout, string fout, double threshold) {

	//0. Prepare genotypes
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes & phenotypes for covariates");
		covariate_engine->residualize(genotype_orig);
		covariate_engine->residualize(phenotype_orig);
	}
	for (int g = 0 ; g < genotype_count ; g ++) genotype_sd[g] = RunningStat(genotype_orig[g]).StandardDeviation();
	for (int p = 0 ; p < phenotype_count ; p ++) phenotype_sd[p] = RunningStat(phenotype_orig[p]).StandardDeviation();
	normalize(genotype_orig);
	normalize(phenotype_orig);

    //0. save the data
	for (int p = 0 ; p < phenotype_count ; p ++) {

		//0.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}

		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

        //0.1 Save to the gene list 
        //0.2 Save the responses for this gene
        //0.3 Save the response names for this gene 

        //0.1 Save the genotype data (X: design matrix) to file 
        string x_out = dout+"/X_" + sutils::int2str(p) + ".npy";
	    if (!futils::createFile(x_out)) LOG.error(x_out + " is impossible to create, check permissions");
		LOG.println("  * Saving eQTL X matirx to [" + x_out + "]");
        const unsigned int x_shape[] = {(int) sample_count,(int) targetGenotypes.size()};
        float* x_data = new float[x_shape[0]*x_shape[1]];
        for(unsigned i=0;i<x_shape[0];i++){
            for(unsigned j=0;j<x_shape[1];j++){
                x_data[i*x_shape[1]+j] = genotype_orig[targetGenotypes[j]][i];
            }
        }
        cnpy::npy_save(x_out,x_data,x_shape,2,"w");
        delete[] x_data;

        //0.2. Save the data to file
        // ofile gdo (gout);
        // gdo << "SAMPLE_ID\t" << phenotype_id[p];
		// for (int g = 0 ; g < targetGenotypes.size() ; g ++)  gdo << "\t" << genotype_id[targetGenotypes[g]];
        // gdo << endl;
        // for (int s = 0 ; s < sample_count; s ++) {
        //     gdo << sample_id[s] << "\t" << sutils::double2str(phenotype_orig[p][s],10);
        //     for (int g = 0 ; g < targetGenotypes.size() ; g ++)    
        //         gdo << "\t" << sutils::double2str(genotype_orig[targetGenotypes[g]][s],10);
        //     gdo << endl;
        // }
        // gdo.close();

		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
    }

	//1. Loop over phenotypes
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		LOG.println("\nProcessing gene [" + phenotype_id[p] + "]");

		//1.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		//1.2. Nominal pass: scan cis-window & compute statistics
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_orig[p]);
			double df = sample_count - 2 - covariate_count;
			double tstat2 = getTstat2(corr, df);
			double pval = getPvalueFromTstat2(tstat2, df);
			double slope = getSlope(corr, phenotype_sd[p], genotype_sd[targetGenotypes[g]]);
			double slope_se = abs(slope) / sqrt(tstat2);
			if (pval <= threshold ) {
				fdo << phenotype_id[p];
				fdo << "\t" << genotype_id[targetGenotypes[g]];
				fdo << "\t" << targetDistances[g];
				fdo << "\t" << pval;
				fdo << "\t" << slope;
				fdo << "\t" << slope_se;
				fdo << endl;
			}
		}
	}
	fdo.close();
}
