#Create input file given parameters

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='Write active noise input file with given input parameters.')
    parser.add_argument('input_dir',
                        help='Directory where input file will be written.')
    #System args
    parser.add_argument('--rho', default='1.414213562373095')
    parser.add_argument('--a', default='1.414213562373095')
    parser.add_argument('--N', default='5')
    parser.add_argument('--La', default='5')
    parser.add_argument('--is_p_x', default='1')
    parser.add_argument('--node_protocol', default='uniform')
    parser.add_argument('--spring_protocol', default='uniform')
    parser.add_argument('--l0', default='1.0')
    parser.add_argument('--K', default='1.0')
    parser.add_argument('--kT', default='0.0')
    parser.add_argument('--potential_type', default='fene')
    parser.add_argument('--drmax', default='0.5')

    #Solver args
    parser.add_argument('--dt', default='0.0001')
    parser.add_argument('--gamma', default='1.0')
    parser.add_argument('--va', default='1.0')
    parser.add_argument('--do_subtract_com', default='0')
    parser.add_argument('--do_pin_node', default='0')

    #Observer args
    parser.add_argument('--output_dir', default='active-network-results')
    parser.add_argument('--network_freq', default='1000')
    parser.add_argument('--thermo_freq', default='1000')

    #LabBench args
    parser.add_argument('--experiment', default='standard')
    parser.add_argument('--equil_steps', default='50000')
    parser.add_argument('--production_steps', default='100000')
    parser.add_argument('--info_freq', default='10000')

    #Active Noise Generator args
    parser.add_argument('--dim', default='1')
    parser.add_argument('--nx', default='32')
    parser.add_argument('--tau', default='1.0')
    parser.add_argument('--Lambda', default='1.0')

    args = parser.parse_args()
    print(args.__dict__)

    try:
        os.makedirs(args.input_dir)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')

    print('Writing input file...')
    com_method = 'com_unconstrained'
    if args.do_subtract_com=='1' and args.do_pin_node=='0':
        com_method = 'com_fixed'
    if args.do_pin_node=='1':
        com_method = 'pinned'

    try:
        os.makedirs(args.input_dir + '/%s' % com_method)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')
    
    with open(args.input_dir + '/%s/%s_net_va=%f_tau=%f_lambda=%f_N=%d_nx=%d.in' % (com_method, args.potential_type, float(args.va),float(args.tau),float(args.Lambda), int(args.N), int(args.nx)), 'w') as f:

        f.write('#System\n')
        f.write('rho = %s\n' % args.rho)
        f.write('a = %s\n' % args.a)
        f.write('N = %s\n' % args.N)
        f.write('La = %s\n' % args.La)
        f.write('is_p_x = %s\n' % args.is_p_x)
        f.write('node_protocol = %s\n' % args.node_protocol)
        f.write('spring_protocol = %s\n' % args.spring_protocol)
        f.write('potential_type = %s\n' % args.potential_type)
        f.write('l0 = %s\n' % args.l0)
        f.write('K = %s\n' % args.K)
        f.write('kT = %s\n' % args.kT)
        f.write('drmax = %s\n' % args.drmax)
        f.write('\n')

        f.write('#Solver\n')
        f.write('dt = %s\n' % args.dt)
        f.write('gamma = %s\n' % args.gamma)
        f.write('va = %s\n' % args.va)
        f.write('do_subtract_com = %s\n' % args.do_subtract_com)
        f.write('do_pin_node = %s\n' % args.do_pin_node)
        f.write('\n')

        f.write('#Observer\n')
        f.write('output_dir = %s\n' % args.output_dir)
        f.write('network_freq = %s\n' % args.network_freq)
        f.write('thermo_freq = %s\n' % args.thermo_freq)
        f.write('\n')

        f.write('#LabBench\n')
        f.write('experiment = %s\n' % args.experiment)
        f.write('equil_steps = %s\n' % args.equil_steps)
        f.write('production_steps = %s\n' % args.production_steps)
        f.write('info_freq = %s\n' % args.info_freq)
        f.write('\n')

        f.write('#Active Noise Generator\n')
        f.write('dim = %s\n' % args.dim)
        f.write('nx = %s\n' % args.nx)
        f.write('tau = %s\n' % args.tau)
        f.write('lambda = %s\n' % args.Lambda)

if __name__=='__main__':
    main()
