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
    string maf_out = dout + "/variants_mafs.txt";
    string g_out = dout + "/gene_names.txt";
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
    string gt_out = dout+"/genotype_data.npy";
    double* gt_data = new double[sample_count*genotype_count];
    const unsigned int gt_shape[] = {sample_count, genotype_count};
    for (unsigned i=0; i<gt_shape[0]; i++){
        for(unsigned j=0; j<gt_shape[1]; j++){
            gt_data[i*gt_shape[1]+j] = (double) genotype_orig[j][i];
        }
    }
    cnpy::npy_save(gt_out,gt_data,gt_shape,2,"w");
    delete[] gt_data;
    LOG.println("Written genotypes to: " + gt_out );
    // Save the original phenotypes to file
    // string pt_out = dout+"/phenotype_data.npy";
    // double* pt_data = new double[sample_count*phenotype_count];
    // const unsigned int pt_shape[] = {sample_count, phenotype_count};
    // for (unsigned i=0; i<pt_shape[0]; i++){
    //     for(unsigned j=0; j<pt_shape[1]; j++){
    //         pt_data[i*pt_shape[1]+j] = (double) phenotype_orig[j][i];
    //     }
    // }
    // cnpy::npy_save(pt_out,pt_data,pt_shape,2,"w");
    // delete[] pt_data;
    // LOG.println("Written genotypes to: " + pt_out );
    // Save indices in genotype for each gene
    ofile g_file (g_out);
    // TODO: use int for v_data
    double* v_data = new double[genotype_count];
    double* y_data = new double[sample_count];
	for (int p = 0 ; p < 1; p ++) {
	// for (int p = 0 ; p < phenotype_count ; p ++) {
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
        if (tg_size == 0) {
            continue; // do not record genes with no variants
        } else {
            g_file << p << "\t" << phenotype_id[p] << "\n"; 
        }
        // Save the responses for this gene
        string v_out = dout+"/v_" + phenotype_id[p] + ".npy";
        string y_out = dout+"/y_" + phenotype_id[p] + ".npy";
        const unsigned int v_shape[] = {tg_size};
        const unsigned int y_shape[] = {sample_count};
        for(unsigned i=0; i<sample_count; i++) y_data[i] = (double) phenotype_orig[p][i];
        for(unsigned j=0; j<tg_size; j++)  v_data[j] = (double) targetGenotypes[j];
        cnpy::npy_save(v_out,v_data,v_shape,1,"w");
        cnpy::npy_save(y_out,y_data,y_shape,1,"w");
		// LOG.println("  * Saved eQTL incides vector to [" + v_out + "]");
		// LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
    }
    delete[] v_data;
    delete[] y_data;
    g_file.close();
    LOG.println("Written genes to: " + g_out); 

    // TODO: remove me after debug
    LOG.error(" DEBUGGING - STOPPED HERE ") ;

    // --------------------------------------------
    

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
        // save the gene name and cis var counts 
        // save the variant names to file
        string v_out = dout+"/Vars_" + phenotype_id[p] + ".txt";
        ofile v_file (v_out);
        for (int g = 0 ; g < targetGenotypes.size() ; g ++)    
            v_file << genotype_id[targetGenotypes[g]] << "\n";
        v_file.close();
    }
    LOG.println("Maximum number of variants in cis = " + sutils::int2str(max_target_genotypes));

    // pass 2: save the data to npy files (used info from pass 1 to allocate memory)
    // LOG.println("\nSaving data (X, y)");
    // double* x_data = new double[sample_count*max_target_genotypes];
    // double* y_data = new double[sample_count];
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
        int tg_size = (int) targetGenotypes.size();
        if (tg_size == 0) continue; // do not record genes with no variants
        //0.2 Save the responses for this gene
        // string y_out = dout+"/y_" + phenotype_id[p] + ".npy";
        // string x_out = dout+"/X_" + phenotype_id[p] + ".npy";
        // const unsigned int x_shape[] = {sample_count,tg_size};
        // const unsigned int y_shape[] = {sample_count};
        // for(unsigned i=0;i<x_shape[0];i++){
            // for(unsigned j=0;j<x_shape[1];j++){
            //     x_data[i*x_shape[1]+j] = (double) genotype_orig[targetGenotypes[j]][i];
            // }
        //     y_data[i] = (double) phenotype_orig[p][i];
        // }
        // cnpy::npy_save(x_out,x_data,x_shape,2,"w");
		// LOG.println("  * Saved eQTL X matirx to [" + x_out + "]");
        // cnpy::npy_save(y_out,y_data,y_shape,1,"w");
		// LOG.println("  * Saved eQTL y vector to [" + y_out + "]");
		// LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
    }
    // delete[] x_data;
    // delete[] y_data;
    // -------------------------------------------------------------------------------------------------------------- 


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

        string p_out = dout + "/pval_" + phenotype_id[p] + ".npy";
        const unsigned int p_shape[] = {tg_size};

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

