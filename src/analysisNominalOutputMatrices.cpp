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

    // Save the data to file 
    // ---------------------
    LOG.println("\nSaving genes and variant information");
    string maf_out = dout + "/meta_variants_mafs.txt";
    string g_out = dout + "/meta_gene_names.txt";
    string dvout = dout + "/cis_vars";
    string dpout = dout + "/cis_pvals";
	if (!futils::createFile(g_out)) LOG.error(g_out + " is impossible to create, check permissions");
	if (!futils::createFile(maf_out)) LOG.error(maf_out + " is impossible to create, check permissions");
    // Save the MAFs to file
    ofile maf_file (maf_out);
    for (int g = 0; g < genotype_count ; g ++) maf_file <<genotype_id[g] << "\t" << genotype_maf[g] << "\n";
    maf_file.close();
    LOG.println("Written mafs to: " + maf_out); 
    // Save the gene names to file
    // Save the original genotypes to file
    // TODO: use int for genotype data
    string gt_f_out = dout+"/genotype_data_float.npy";
    double* gt_f_data = new double[sample_count*genotype_count];
    const unsigned int gt_f_shape[] = {(unsigned int) sample_count, (unsigned int) genotype_count};
    for (unsigned i=0; i<gt_f_shape[0]; i++){
        for(unsigned j=0; j<gt_f_shape[1]; j++){
            gt_f_data[i*gt_f_shape[1]+j] = (double) genotype_orig[j][i];
        }
    }
    cnpy::npy_save(gt_f_out,gt_f_data,gt_f_shape,2,"w");
    delete[] gt_f_data;
    LOG.println("Written genotypes to: " + gt_f_out );

    // int version 
    string gt_out = dout+"/genotype_data_int.npy";
    int64_t* gt_data = new int64_t[sample_count*genotype_count];
    const unsigned int gt_shape[] = {(unsigned int) sample_count, (unsigned int) genotype_count};
    int64_t gt_value; // should be 0, 1, 2
    for (unsigned i=0; i<gt_shape[0]; i++){
        for(unsigned j=0; j<gt_shape[1]; j++){
            // because there are imputed genotype values, do some rounding here
            gt_value = (int64_t) round(genotype_orig[j][i]);
            gt_value = (gt_value > 2) ? 2 : gt_value;
            gt_value = (gt_value < 0) ? 0 : gt_value;
            gt_data[i*gt_shape[1]+j] = gt_value;
        }
    }
    cnpy::npy_save(gt_out,gt_data,gt_shape,2,"w");
    delete[] gt_data;
    LOG.println("Written genotypes to: " + gt_out);

    // Save indices in genotype for each gene
    ofile g_file (g_out);
    int64_t* v_data = new int64_t[genotype_count];
    int* used_phenotype_index = new int[phenotype_count];
    int used_phenotype_count = 0;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		// Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
        int tg_size = (int) targetGenotypes.size();
        if (tg_size == 0) { continue; }// do not record genes with no variants 
        g_file << phenotype_id[p] << "\n"; 
        // Save the responses for this gene
        string v_out = dvout+"/v_" + phenotype_id[p] + ".npy";
        const unsigned int v_shape[] = {(unsigned int) tg_size};
        for(unsigned j=0; j<tg_size; j++)  v_data[j] = (int64_t) targetGenotypes[j];
        cnpy::npy_save(v_out,v_data,v_shape,1,"w");
        used_phenotype_index[used_phenotype_count] = p;
        used_phenotype_count++;
    }
    LOG.println("Written genes to: " + g_out); 

    string pt_out = dout+"/phenotype_data.npy";
    double* pt_data = new double[sample_count*used_phenotype_count];
    const unsigned int pt_shape[] = {(unsigned int) sample_count, (unsigned int) used_phenotype_count};
    for (unsigned i=0; i<pt_shape[0]; i++){
        for(unsigned j=0; j<pt_shape[1]; j++){
            pt_data[i*pt_shape[1]+j] = (double) phenotype_orig[used_phenotype_index[j]][i];
        }
    }
    cnpy::npy_save(pt_out,pt_data,pt_shape,2,"w");
    delete[] pt_data;
    delete[] used_phenotype_index;
    delete[] v_data;
    g_file.close();

    // --------------------------------------------
	//0. Prepare genotypes
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes & phenotypes for covariates");
		covariate_engine->residualize(genotype_orig);
		covariate_engine->residualize(phenotype_orig);
        // print covariates to file
        // TODO: maybe remve later
        string cov_out = dout+"/covariate_data_float.npy";
        double* cov_data = new double[sample_count*(covariate_count+1)];
        const unsigned int cov_shape[] = {(unsigned int) sample_count, (unsigned int) covariate_count+1};
        for (unsigned i=0; i<cov_shape[0]; i++){
            for(unsigned j=0; j<cov_shape[1]; j++){
                cov_data[i*cov_shape[1]+j] = (double) covariate_engine->covarM(i, j);
            }
        }
        cnpy::npy_save(cov_out,cov_data,cov_shape,2,"w");
        delete[] cov_data;
        LOG.println("Written covariates to: " + cov_out );
	}
	for (int g = 0 ; g < genotype_count ; g ++) genotype_sd[g] = RunningStat(genotype_orig[g]).StandardDeviation();
	for (int p = 0 ; p < phenotype_count ; p ++) phenotype_sd[p] = RunningStat(phenotype_orig[p]).StandardDeviation();

    // *** save regressed genotypes for debugging TODO: maybe remove later
    string ro_out = dout+"/genotype_data_reg_out_float.npy";
    double* ro_data = new double[sample_count*genotype_count];
    const unsigned int ro_shape[] = {(unsigned int) sample_count, (unsigned int) genotype_count};
    for (unsigned i=0; i<ro_shape[0]; i++){
        for(unsigned j=0; j<ro_shape[1]; j++){
            ro_data[i*ro_shape[1]+j] = (double) genotype_orig[j][i];
        }
    }
    cnpy::npy_save(ro_out,ro_data,ro_shape,2,"w");
    delete[] ro_data;
    LOG.println("Written regressed genotypes to: " + ro_out );

    // TODO: remove me after debug
    // LOG.error(" DEBUGGING - STOPPED HERE ") ;
    // -------------------------------------------------------------------------------------------------------------- 
    // pass 1: count how many variants there are per gene, save the variants to file as well

    int max_target_genotypes = 0;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		// enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}

        int tg_size = (int) targetGenotypes.size();
		LOG.println("  * Gene: " + phenotype_id[p] + "\t Number of variants in cis = " + sutils::int2str(tg_size));
        if (tg_size == 0) continue; // do not record genes with no variants

        max_target_genotypes = max(max_target_genotypes, tg_size);
    }
    LOG.println("Maximum number of variants in cis = " + sutils::int2str(max_target_genotypes));

	normalize(genotype_orig);
	normalize(phenotype_orig);

    double* p_data = new double[max_target_genotypes];
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
        int tg_size = (int) targetGenotypes.size();
        if (tg_size == 0) continue; // do not record genes with no variants

        string p_out = dpout + "/pval_" + phenotype_id[p] + ".npy";
        const unsigned int p_shape[] = {(unsigned int) tg_size};

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
                p_data[g] = pval;
			}
		}
        cnpy::npy_save(p_out,p_data,p_shape,1,"w");
		LOG.println("  * Saved eQTL p-value vector to [" + p_out + "]");
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
    delete[] p_data;
	fdo.close();
}

