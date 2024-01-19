#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Intersect readID-enhancerBC-plasmidBC file with enhancerID-enhancerBC to return enhancerID-plasmidBC assignments.",)
    parser.add_argument('--enhancerBC', '-ebc', type=str, required=True,
                        help='Path to enhancerID-enhancerBC file (a table separated by tab)')
    parser.add_argument('--readBC', '-rbc', type=str, required=True,
                        help='Path to readID-enhancerBC-plasmidBC file (a table separated by tab)')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Path to output file (a table with enhancerID-plasmidBC assignments separated by tab).')
    parser.add_argument('--ratio_threshold', '-rthr', type=int, required=False, default = 10,
                        help='Ratio threshold')
    parser.add_argument('--contribution_threshold', '-cthr', type=float, required=False, default = 0.95,
                        help='Contribution threshold')
    return parser

def main():
    """
    The main executable function
    """
    parser = make_argument_parser()
    args = parser.parse_args()
    ebc = args.enhancerBC
    rbc = args.readBC
    rthr = args.ratio_threshold
    cthr = args.contribution_threshold
    o = args.output
    
    id2bce=pd.read_csv(ebc, sep='\t', header=None, names=['enhancerID', 'enhancerBarcode'])
    bcs = pd.read_csv(rbc, sep='\t', header=None, names=['readID', 'plasmidBarcode', 'enhancerBarcode'])
    id2csbc = pd.merge(id2bce, bcs, on='enhancerBarcode')
    id2csbc_nd = id2csbc[['enhancerID', 'plasmidBarcode']].drop_duplicates()
    print('Number of unique assignments:', id2csbc_nd.shape[0])
    id2csbc_un = id2csbc_nd.drop_duplicates(subset=['plasmidBarcode'], keep = False)
    id2csbc_all = id2csbc_nd.drop_duplicates(subset=['plasmidBarcode'])
    print('Number of unique plasmid barcodes:', id2csbc_all.shape[0])
    print('Number of unique plasmid barcodes uniquely assigned to an enhancer:', id2csbc_un.shape[0])
    print('Complexity:', id2csbc_un.shape[0]/id2csbc_nd.shape[0])
    id2csbc_un.to_csv(o, header=None, sep='\t', index=False)
    id2csbc_d = id2csbc[~id2csbc['plasmidBarcode'].isin(id2csbc_un['plasmidBarcode'])]
    id2csbc_d = id2csbc_d[['enhancerID', 'plasmidBarcode']]
    id2csbc_dg = id2csbc_d.groupby(["enhancerID", "plasmidBarcode"]).size().reset_index()
    id2csbc_dg.columns = ["enhancerID","plasmidBarcode", 'counts']
    id2csbc_dg['Combined_name'] = id2csbc_dg['enhancerID'] + '_' + id2csbc_dg['plasmidBarcode']
    most_common = id2csbc_dg.sort_values(['counts'], ascending=False).drop_duplicates(["plasmidBarcode"], keep='first')
    second_most_common = id2csbc_dg[~id2csbc_dg['Combined_name'].isin(most_common['Combined_name'])].sort_values(['counts'], ascending=False).drop_duplicates(["plasmidBarcode"], keep='first')
    id2csbc_dg = most_common.merge(second_most_common, on='plasmidBarcode')  
    id2csbc_dg['Ratio'] = id2csbc_dg['counts_x']/id2csbc_dg['counts_y'] 
    id2csbc_dg = id2csbc_dg[['enhancerID_x', 'counts_x', 'plasmidBarcode', 'Ratio']]
    id2csbc_dg.columns = ['enhancerID', 'counts', 'plasmidBarcode', 'Ratio']
    total_counts = id2csbc_d.groupby(["plasmidBarcode"]).size().reset_index()
    total_counts.columns = ['plasmidBarcode', 'total_counts']
    id2csbc_dg = id2csbc_dg.merge(total_counts, on='plasmidBarcode')
    id2csbc_dg['Contribution'] = id2csbc_dg['counts']/id2csbc_dg['total_counts']
    id2csbc_dg = id2csbc_dg[['enhancerID', 'plasmidBarcode', 'Ratio', 'Contribution']]
    print('Number of uniquely assigned plasmid barcodes recovered with contribution > ', cthr, ' : ', sum(id2csbc_dg['Contribution'] > cthr))
    sub = id2csbc_dg[id2csbc_dg['Contribution'] > cthr]
    c = id2csbc_d[~id2csbc_d['plasmidBarcode'].isin(sub['plasmidBarcode'])].drop_duplicates()
    total = pd.concat([c, id2csbc_un, sub[['enhancerID', 'plasmidBarcode']]])
    num = pd.concat([id2csbc_un, sub[['enhancerID', 'plasmidBarcode']]])
    print('Number of unique assignments after filtering:', total.shape[0])
    print('Number of unique plasmid barcodes uniquely assigned to an enhancer after filtering:', num.shape[0])
    print('Complexity after filtering:', num.shape[0]/total.shape[0])
    print('Number of uniquely assigned plasmid barcodes recovered with ratio > ',  rthr, ' : ', sum(sub['Ratio'] > rthr))
    sub = sub[sub['Ratio'] > rthr]
    c = id2csbc_d[~id2csbc_d['plasmidBarcode'].isin(sub['plasmidBarcode'])].drop_duplicates()
    total = pd.concat([c, id2csbc_un, sub[['enhancerID', 'plasmidBarcode']]])
    num = pd.concat([id2csbc_un, sub[['enhancerID', 'plasmidBarcode']]])
    print('Number of unique assignments after filtering:', total.shape[0])
    print('Number of unique plasmid barcodes uniquely assigned to an enhancer after filtering:', num.shape[0])
    print('Complexity after filtering:',  num.shape[0]/total.shape[0])
    id2csbc_dg.to_csv(o.rsplit('/', 1)[0] + '/' + 'Recovered_stats_' + o.rsplit('/', 1)[1], header=None, sep='\t', index=False)
    num.to_csv(o.rsplit('/', 1)[0] + '/' + 'Recovered_' + o.rsplit('/', 1)[1], header=None, sep='\t', index=False)
     
if __name__ == "__main__":
    main()
