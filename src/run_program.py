import os, subprocess, sys
sys.path.append('/global_utils/src/')
# sys.path.append('global_utils/src/')
import module_utils
import goqc
        
def run_program( arg_list ):
    """
    Parameters:
    -i <input_go_files> - can take a list of files
    -type <input_go_file_type> - default davidgo
    -name <analysis_name> - default go
    -o <output_dir>
    -pvalue <pvalue_cutoff> - pvalue cutoff for showing a term on GO plot - default 0.3
    """
    print('ARG LIST: {}'.format(str(arg_list)))
    input_args = module_utils.getArgument( arg_list, '-i', 'list' )
    output_dir = module_utils.getArgument( arg_list, '-o' )
    go_file_type = module_utils.getArgument( arg_list, '-type', 'implicit', 'davidgo' )
    go_name = module_utils.getArgument( arg_list, '-name', 'implicit', 'go' )
    pvalue = float(module_utils.getArgument( arg_list, '-pvalue', 'implicit', 0.3 ))
    if output_dir not in [[], '']:
        os.chdir( output_dir )
    if input_args != []:
        # create input JSON
        input_json = {'analysis_name': go_name, 'input_file': input_args, 'input_file_type': go_file_type, \
                      'output_dir': output_dir, 'pvalue_cutoff': pvalue}
        # run GO QC
        goqc.goqc( input_json )
    return

if __name__ == '__main__':
    print('in run_program.py')
    run_program( sys.argv[1:] )
