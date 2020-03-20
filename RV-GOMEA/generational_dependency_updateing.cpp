//
// Created by chantal on 31-3-19.
//


double **dependency_matrix                = (double **) Malloc( number_of_parameters*sizeof( double * ) );
double **checked_pairs                = (double **) Malloc( number_of_parameters*sizeof( double * ) );
double *first_individual                 = (double *) Malloc( number_of_parameters*sizeof( double ) );
double *second_individual                = (double *) Malloc( number_of_parameters*sizeof( double ) );
double *fitness_of_first_individual      = (double *) Malloc( (number_of_parameters + 1)*sizeof( double * ) );
double *dependency_pairs = (int **) Malloc( ((number_of_parameters*number_of_parameters)/2)*sizeof(int * ) );
int number_of_parameters = l;
int current_waiting_position;
int sparse_learning;



void initialize(){
    // rv gomea flags
    static_linkage_tree = 1; dependency_learning = 1; evolve_learning = 1; pruning_ub = 100; continued_learning=1;
    // to create a marginal product FOS we need to have a sparse tree, this flag should be one
    int sparse_tree = 1;
    for(int j = 0; j < number_of_parameters; j++ ){
        dependency_matrix[j] = (double *) Malloc( number_of_parameters*sizeof( double ) );
        checked_matrix[j] = (int *) Malloc( number_of_parameters*sizeof( int ) );
    }
    pairs_per_run = number_of_parameters;
    number_of_checked_pairs = 0;
    int counter = 0;
    for (i = 0; i < number_of_parameters; i++) {
        for (j = i + 1; j < number_of_parameters; j++) {
            // add pairs to evaluate to the list
            dependency_pairs[counter] = (int *) Malloc(2 * sizeof(int));
            dependency_pairs[counter][0] = i;
            dependency_pairs[counter][1] = j;
            counter++;
        }
    }
    number_of_pairs = counter;
    for (int i = counter - 1; i >= 0; --i) {
        //generate a random number [0, n-1]
        int j = randomInt(i+1);

        //swap the last element with element at random index
        int *temp = dependency_pairs[i];
        dependency_pairs[i] = dependency_pairs[j];
        dependency_pairs[j] = temp;
    }

    // fill matrix with 0s
    for (i = 0; i < number_of_parameters; i++) {
        for (j = i; j < number_of_parameters; j++) {
            dependency_matrix[i][j] = 0.0;
            dependency_matrix[j][i] = 0.0;
            checked_matrix[i][j] = 0;
            checked_matrix[j][i] = 0;
        }
    }
    current_waiting_position = 0;
}

FOS *learnFOS(population_index){
    if( current_waiting_position == 0 ){
        if( current_waiting_position == 0 ){
            if(number_of_generations[population_index] != 0){
                ezilaitiniCovarianceMatrices(population_index);
                ezilaitiniFOS(linkage_model[population_index]);
            }
            linkage_model[population_index] = learnLinkageTreeRVGOMEA( population_index );
            initializeCovarianceMatrices( population_index );

        }

        if( number_of_generations[population_index] == 0 ) {
            initializeDistributionMultipliers( population_index );
        }

    }
    else if (current_waiting_position > 0){
        current_waiting_position -= 1;
    }
}


FOS *learnLinkageTreeRVGOMEA( int population_index )
{
    int i;
    FOS *new_FOS;
    if( evolve_learning ){
        evolveDifferentialDependencies( population_index );
    }
    new_FOS = learnLinkageTree( full_covariance_matrix[population_index], dependency_matrix, checked_matrix);
    return( new_FOS );
}

void evolveDifferentialDependencies( int population_index ) {
    int i, j, k;
    double *individual_to_compare = (double *) Malloc(number_of_parameters * sizeof(double));
    double constraint_value;

    // initialize if no pairs are checked yet
    if (number_of_checked_pairs == 0) {
        double rand = randomRealUniform01();
        rand = 0.7;

        for (k = 0; k < number_of_parameters; k++) {
            double min = lower_init_ranges[k], max = upper_init_ranges[k];
            getMinMaxofPopulation(k, population_index, &min, &max);
            if (nround(min, 2) == nround(max, 2)) {
                max = upper_init_ranges[k];
            }
            first_individual[k] = min + ((max - min) * rand * 0.5);
            double parameter_diff = (max - min) * 0.5 * rand;
            second_individual[k] = parameter_diff + first_individual[k];
            individual_to_compare[k] = first_individual[k];
        }

        double objective_value, old_constraint, old_objective;
        // fill evaluation storage
        installedProblemEvaluation(problem_index, first_individual, &(old_objective), &(old_constraint),
                                   number_of_parameters, NULL, NULL, 0, 0);
        differential_grouping_evals = 1+ number_of_parameters;
        fitness_of_first_individual[number_of_parameters] = old_objective;
        fitness_of_first_individual[0] = old_objective;
        for (k = 0; k < number_of_parameters; k++) {
            individual_to_compare[k] = second_individual[k];
            installedProblemEvaluation(problem_index, individual_to_compare, &(objective_value), &(constraint_value), 1, &(k), &(first_individual[k]), old_objective, old_constraint);

            fitness_of_first_individual[k] = objective_value;
            individual_to_compare[k] = first_individual[k];
        }
        int counter = number_of_pairs;
        for (int i = counter - 1; i >= 0; --i) {
            int j = randomInt(i+1);

            //swap the last element with element at random index
            int *temp = dependency_pairs[i];
            dependency_pairs[i] = dependency_pairs[j];
            dependency_pairs[j] = temp;
        }

    } else {
        for (k = 0; k < number_of_parameters; k++) {
            individual_to_compare[k] = first_individual[k];
        }
    }

    iteration += 1;
    int max_index = number_of_checked_pairs + pairs_per_run;
    if (max_index >= number_of_pairs) {
        max_index = number_of_pairs;
    }

    double original_objective = fitness_of_first_individual[number_of_parameters];

    for (k = 0; k < number_of_parameters; k++) {
        individual_to_compare[k] = first_individual[k];
    }
    int found_dependencies = 0;
    double max_dependency = 0.0;
    for (k = number_of_checked_pairs; k < max_index; k++) {
        i = dependency_pairs[k][0];
        j = dependency_pairs[k][1];

        double change_i, change_j, change_i_j;
        change_i = fitness_of_first_individual[i];
        change_j = fitness_of_first_individual[j];

        individual_to_compare[i] = second_individual[i];
        individual_to_compare[j] = second_individual[j];
        installedProblemEvaluation(temp_problem_index, individual_to_compare, &(change_i_j), &(constraint_value),
                                   1, &(j), &(first_individual[j]), fitness_of_first_individual[i], 0);
        differential_grouping_evals+=1;
        individual_to_compare[i] = first_individual[i];
        individual_to_compare[j] = first_individual[j];

        double delta_i, delta_j;

        change_i = change_i/original_objective;
        change_j = change_j/original_objective;


        delta_i = fabs(1.0 - change_i);
        delta_j = fabs(change_j - change_i_j);


        delta_i = nround(delta_i, 12);
        delta_j = nround(delta_j, 12);

        double dependency = 0.0;
        double inverted_difference;

        if(delta_j == 0.0) {
            double temp = delta_i;
            delta_i = delta_j;
            delta_j = temp;
        }
        if(delta_j != 0.0){
            inverted_difference = fabs(delta_i/delta_j);
            if(inverted_difference > 1.0){
                inverted_difference = fabs((double)delta_j/delta_i);
            }
        } else{
            inverted_difference = 1.0;
        }
        dependency = 1-inverted_difference;
        if (inverted_difference < 1) {
            found_dependencies += 1;
        } else{
            dependency = 0.0;
        }
        dependency_matrix[i][j] = dependency;
        dependency_matrix[j][i] = dependency;

        max_dependency = max(dependency, max_dependency);
        checked_matrix[i][j] = 1;
        checked_matrix[j][i] = 1;
    }
    total_dependencies_found += found_dependencies;
    number_of_checked_pairs += pairs_per_run;
    if (found_dependencies == 0) {
        int found_dependencies_per_run = total_dependencies_found / iteration;
        if (found_dependencies_per_run < minimal_dependencies_per_run) {
            current_waiting_position = number_of_waiting_cycles;
            number_of_waiting_cycles *= 2;
            iteration = 0; total_dependencies_found = 0;
        }
    }
    if (number_of_checked_pairs >= number_of_pairs){
        number_of_checked_pairs = 0;
        current_waiting_position = number_of_waiting_cycles;
        number_of_waiting_cycles *= 2;
        iteration = 0; total_dependencies_found = 0;
    }

    free(individual_to_compare);
}


FOS *learnLinkageTree( double **covariance_matrix , double **dependency_matrix, int **checked_matrix)
{
    char     done;
    int      i, j, r0, r1, rswap, *indices, *order, *sorted,
            FOS_index, **mpm, mpm_length,
            **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
            *NN_chain, NN_chain_length;
    double   mul0, mul1, **MI_matrix;
    FOS *new_FOS;

    /* Compute Mutual Information matrix */
    MI_matrix = NULL;
    if( learn_linkage_tree && !dependency_learning)
        MI_matrix = computeMIMatrix( covariance_matrix, number_of_parameters );

    /* Initialize MPM to the univariate factorization */
    order                 = randomPermutation( number_of_parameters );
    mpm                   = (int **) Malloc( number_of_parameters*sizeof( int * ) );
    mpm_number_of_indices = (int *) Malloc( number_of_parameters*sizeof( int ) );
    mpm_length            = number_of_parameters;
    mpm_new               = NULL;
    for( i = 0; i < number_of_parameters; i++ )
    {
        indices                  = (int *) Malloc( 1*sizeof( int ) );
        indices[0]               = order[i];
        mpm[i]                   = indices;
        mpm_number_of_indices[i] = 1;
    }
    free( order );

    /* Initialize LT to the initial MPM */
    new_FOS                     = (FOS*) Malloc(sizeof(FOS));
    new_FOS->length             = number_of_parameters+number_of_parameters-1;
    new_FOS->sets               = (int **) Malloc( new_FOS->length*sizeof( int * ) );
    new_FOS->set_length         = (int *) Malloc( new_FOS->length*sizeof( int ) );
    FOS_index                                   = 0;
    for( i = 0; i < mpm_length; i++ )
    {
        new_FOS->sets[FOS_index]       = mpm[i];
        new_FOS->set_length[FOS_index] = mpm_number_of_indices[i];
        FOS_index++;
    }

    /* Initialize similarity matrix */
    S_matrix = NULL;
    if( !random_linkage_tree ){
        S_matrix = (double **) Malloc( number_of_parameters*sizeof( double * ) );
        for( i = 0; i < number_of_parameters; i++ )
            S_matrix[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );
    }

    if( learn_linkage_tree )
    {
        if ( dependency_learning ) {
            for (i = 0; i < mpm_length; i++)
                for (j = 0; j < mpm_length; j++)
                    S_matrix[i][j] = dependency_matrix[mpm[i][0]][mpm[j][0]];
            for (i = 0; i < mpm_length; i++)
                S_matrix[i][i] = 0;
        }
        else{
            for( i = 0; i < mpm_length; i++ )
                for( j = 0; j < mpm_length; j++ )
                    S_matrix[i][j] = MI_matrix[mpm[i][0]][mpm[j][0]];
            for( i = 0; i < mpm_length; i++ )
                S_matrix[i][i] = 0;

            for( i = 0; i < number_of_parameters; i++ )
                free( MI_matrix[i] );
            free( MI_matrix );
        }

    }
    else if( random_linkage_tree )
    {
        S_vector = (double *) Malloc( number_of_parameters*sizeof(double));
        for( i = 0; i < number_of_parameters; i++ )
            S_vector[i] = randomRealUniform01();
    }
    else if( static_linkage_tree )
    {
        if( problem_index == 105 || problem_index == 106 )
        {
            for( i = 0; i < number_of_parameters-1; i++ )
            {
                for( j = i+1; j < number_of_parameters; j++ )
                {
                    S_matrix[i][j] = 1.0 / covariance_matrix[mpm[i][0]][mpm[j][0]];
                    S_matrix[j][i] = S_matrix[i][j];
                }
                S_matrix[i][i] = 0.0;
            }
        }
        else if (dependency_learning){
            for (i = 0; i < mpm_length; i++)
                for (j = 0; j < mpm_length; j++)
                    S_matrix[i][j] = dependency_matrix[mpm[i][0]][mpm[j][0]];
            for (i = 0; i < mpm_length; i++)
                S_matrix[i][i] = 0;
        }
        else
        {
            for( i = 0; i < mpm_length; i++ )
            {
                for( j = 0; j < i; j++ )
                {
                    if( mpm[i][0] < block_start || mpm[j][0] < block_start ) S_matrix[i][j] = randomRealUniform01();
                    else if( (mpm[i][0]-block_start)/block_size == (mpm[j][0]-block_start)/block_size ) S_matrix[i][j] = randomRealUniform01() + 1e8;
                    else S_matrix[i][j] = randomRealUniform01() + 1e3;
                    S_matrix[j][i] = S_matrix[i][j];
                }
                S_matrix[i][i] = 0;
            }
        }
    }

    int *keep_FOS_element = (int *) Malloc( ((number_of_parameters*2))*sizeof( int ) );
    for( i = 0; i < number_of_parameters*2; i++ ){
        keep_FOS_element[i] = 1;
    }

    NN_chain        = (int *) Malloc( (number_of_parameters+2)*sizeof( int ) );
    NN_chain_length = 0;
    done            = 0;
    while( !done )
    {
        if( NN_chain_length == 0 )
        {
            NN_chain[NN_chain_length] = randomInt( mpm_length );
            NN_chain_length++;
        }

        if( NN_chain[NN_chain_length-1] >= mpm_length ) NN_chain[NN_chain_length-1] = mpm_length-1;

        while( NN_chain_length < 3 )
        {
            NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
            NN_chain_length++;
        }

        while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
        {
            NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
            if( ((getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length]) == getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length-2])))
                && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
                NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
            NN_chain_length++;
            if( NN_chain_length > number_of_parameters ){
                break;
            }
        }
        r0 = NN_chain[NN_chain_length-2];
        r1 = NN_chain[NN_chain_length-1];

        if( r1 >= mpm_length || r0 >= mpm_length || mpm_number_of_indices[r0]+mpm_number_of_indices[r1] > FOS_element_ub )
        {
            NN_chain_length = 1;
            NN_chain[0] = 0;
            if( FOS_element_ub < number_of_parameters )
            {
                done = 1;
                for( i = 1; i < mpm_length; i++ )
                {
                    if( mpm_number_of_indices[i] + mpm_number_of_indices[NN_chain[0]] <= FOS_element_ub ) done = 0;
                    if( mpm_number_of_indices[i] < mpm_number_of_indices[NN_chain[0]] ) NN_chain[0] = i;
                }
                if( done ) break;
            }
            continue;
        }

        if( r0 > r1 )
        {
            rswap = r0;
            r0    = r1;
            r1    = rswap;
        }
        NN_chain_length -= 3;

        if( r1 < mpm_length && r1 != r0 ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
        {
            indices = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );

            i = 0;
            for( j = 0; j < mpm_number_of_indices[r0]; j++ )
            {
                indices[i] = mpm[r0][j];
                i++;
            }
            for( j = 0; j < mpm_number_of_indices[r1]; j++ )
            {
                indices[i] = mpm[r1][j];
                i++;
            }

            new_FOS->sets[FOS_index] = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );
            new_FOS->set_length[FOS_index] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
            sorted = mergeSortInt(indices, mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
            for( j = 0; j < mpm_number_of_indices[r0]+mpm_number_of_indices[r1]; j++ ){
                new_FOS->sets[FOS_index][j] = indices[sorted[j]];
            }

            if(getSimilarity(r0, r1) <= epsilon){
                keep_FOS_element[FOS_index] = 0;
            }

            if( mpm_number_of_indices[r0]+mpm_number_of_indices[r1] > pruning_ub && keep_FOS_element[FOS_index] ){
                keep_FOS_element[FOS_index] = 0;
            }
            if( keep_FOS_element[FOS_index] ){
                // we know we will merge r0 and r1, now lets check if they are all completely dependent
                int completely_dependent = 1;
                int all_checked = 1;
                for (i = 0; i < mpm_number_of_indices[r0]; i++){
                    for (j = 0; j< mpm_number_of_indices[r1]; j++){
                        if (dependency_matrix[mpm[r0][i]][mpm[r1][j]] <= 0.0){
                            if(checked_matrix[mpm[r0][i]][mpm[r1][j]]==0){
                                all_checked = 0;
                            }
                            completely_dependent = 0;
                            break;
                        }
                    }
                }
                if ( completely_dependent ) { // remove subsets that build this set
                    //remove r1
                    int first_set_element = mpm[r0][0];
                    int set_length = mpm_number_of_indices[r0];
                    for (i = 0; i < FOS_index; i++) {
                        if (new_FOS->set_length[i] == set_length && new_FOS->sets[i][0] == first_set_element) {
                            keep_FOS_element[i] = 0;
                        }
                    }
                    //remove r0
                    first_set_element = mpm[r1][0];
                    set_length = mpm_number_of_indices[r1];
                    for (i = 0; i < FOS_index; i++) {
                        if (new_FOS->set_length[i] == set_length && new_FOS->sets[i][0] == first_set_element) {
                            keep_FOS_element[i] = 0;
                        }
                    }
                }
                if( all_checked ){
                    keep_FOS_element[FOS_index] = 0;
                }

            }


            free( sorted );
            free( indices );

            mul0 = ((double) mpm_number_of_indices[r0])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
            mul1 = ((double) mpm_number_of_indices[r1])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
            if( random_linkage_tree )
            {
                S_vector[r0] = mul0*S_vector[r0]+mul1*S_vector[r1];
            }
            else
            {
                for( i = 0; i < mpm_length; i++ )
                {
                    if( (i != r0) && (i != r1) )
                    {
                        S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
                        S_matrix[r0][i] = S_matrix[i][r0];
                    }
                }
            }

            mpm_new                   = (int **) Malloc( (mpm_length-1)*sizeof( int * ) );
            mpm_new_number_of_indices = (int *) Malloc( (mpm_length-1)*sizeof( int ) );
            mpm_new_length            = mpm_length-1;
            for( i = 0; i < mpm_new_length; i++ )
            {
                mpm_new[i]                   = mpm[i];
                mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
            }

            mpm_new[r0]                   = new_FOS->sets[FOS_index];
            mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
            if( r1 < mpm_length-1 )
            {
                mpm_new[r1]                   = mpm[mpm_length-1];
                mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length-1];

                if( random_linkage_tree )
                {
                    S_vector[r1] = S_vector[mpm_length-1];
                }
                else
                {
                    for( i = 0; i < r1; i++ )
                    {
                        S_matrix[i][r1] = S_matrix[i][mpm_length-1];
                        S_matrix[r1][i] = S_matrix[i][r1];
                    }

                    for( j = r1+1; j < mpm_new_length; j++ )
                    {
                        S_matrix[r1][j] = S_matrix[j][mpm_length-1];
                        S_matrix[j][r1] = S_matrix[r1][j];
                    }
                }
            }

            for( i = 0; i < NN_chain_length; i++ )
            {
                if( NN_chain[i] == mpm_length-1 )
                {
                    NN_chain[i] = r1;
                    break;
                }
            }

            free( mpm );
            free( mpm_number_of_indices );
            mpm                   = mpm_new;
            mpm_number_of_indices = mpm_new_number_of_indices;
            mpm_length            = mpm_new_length;

            if( mpm_length == 1 ){
                done = 1;
            }

            FOS_index++;
        }
    }

    i = 0;
    j = 0;
    int new_lenght = 0;
    while (i < FOS_index){
        if (keep_FOS_element[i]) {
            if (i > j){
                while(keep_FOS_element[j] && j < FOS_index )
                    j ++;
                new_FOS->set_length[j] = new_FOS->set_length[i];
                free(new_FOS->sets[j]);
                new_FOS->sets[j] = (int *) Malloc( (new_FOS->set_length[j])*sizeof( int ) );
                for(int k = 0; k < new_FOS->set_length[j]; k++ ){
                    new_FOS->sets[j][k] = new_FOS->sets[i][k];
                }
                keep_FOS_element[j] = 1;
                keep_FOS_element[i] = 0;
            }
            if( i == j ){
                j++;
            }
            new_lenght += 1;
        }
        i ++;
    }
    for(i = new_lenght; i<FOS_index; i++){
        free(new_FOS->sets[i]);
    }
    FOS_index = new_lenght;


    new_FOS->length = FOS_index;

    free( NN_chain );

    free( mpm_new );
    free( mpm_number_of_indices );
    free( keep_FOS_element );

    if( random_linkage_tree )
        free( S_vector );
    else
    {
        for( i = 0; i < number_of_parameters; i++ )
            free( S_matrix[i] );
        free( S_matrix );
    }

    return( new_FOS );
}

void inheritDistributionMultipliers( FOS *new_FOS, FOS *prev_FOS, double *multipliers )
{
    int      i, *permutation;
    double   *multipliers_copy;

    multipliers_copy = (double*) Malloc(new_FOS->length*sizeof(double));
    for( i = 0; i < new_FOS->length; i++ )
        multipliers_copy[i] = multipliers[i];
    permutation = matchFOSElements( new_FOS, prev_FOS );

    for( i = 0; i < number_of_parameters; i++ )
        multipliers[permutation[i]] = multipliers_copy[i];


    for( i = 0; i < new_FOS->length; i++ ){
        multipliers[i] = multipliers_copy[permutation[i]];
    }

    free( multipliers_copy );
    free( permutation );
}


int *matchFOSElements( FOS *new_FOS, FOS *prev_FOS )
{
    int      i, j, a, b, matches, *permutation, *hungarian_permutation,
            **FOS_element_similarity_matrix;

    permutation = (int *) Malloc( new_FOS->length*sizeof(int));
    FOS_element_similarity_matrix = (int**) Malloc((prev_FOS->length-number_of_parameters)*sizeof(int*));
    for( i = 0; i < prev_FOS->length-number_of_parameters; i++ )
        FOS_element_similarity_matrix[i] = (int*) Malloc((new_FOS->length-number_of_parameters)*sizeof(int));

    for( i = 0; i < number_of_parameters; i++ )
    {
        for( j = 0; j < number_of_parameters; j++ )
        {
            if( prev_FOS->sets[i][0] == new_FOS->sets[j][0] )
            {
                permutation[i] = j;
                break;
            }
        }
    }
    for( i = number_of_parameters; i < prev_FOS->length; i++ )
    {
        for( j = number_of_parameters; j < new_FOS->length; j++ )
        {
            a = 0; b = 0;
            matches = 0;
            while( a < prev_FOS->set_length[i] && b < new_FOS->set_length[j] )
            {
                if( prev_FOS->sets[i][a] < new_FOS->sets[j][b] )
                {
                    a++;
                }
                else if( prev_FOS->sets[i][a] > new_FOS->sets[j][b] )
                {
                    b++;
                }
                else
                {
                    a++;
                    b++;
                    matches++;
                }
            }
            FOS_element_similarity_matrix[i-number_of_parameters][j-number_of_parameters] = (int) 10000*(2.0*matches/(prev_FOS->set_length[i]+new_FOS->set_length[j]));
        }
    }

    for( i = 0; i < new_FOS->length; i++ )
    {
        int max_index = 0;
        int max_similarity = -1;
        for( j = number_of_parameters; j < prev_FOS->length; j++ )
        {
            if(FOS_element_similarity_matrix[j-number_of_parameters][i-number_of_parameters]>max_similarity){
                max_index = j;
                max_similarity = FOS_element_similarity_matrix[j-number_of_parameters][i-number_of_parameters];
            }
        }
        permutation[i] = max_index;
    }

    for( i = 0; i < prev_FOS->length-number_of_parameters; i++ )
        free( FOS_element_similarity_matrix[i] );
    free( FOS_element_similarity_matrix );

    return( permutation );
}

