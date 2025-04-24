#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2025
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2025-03-02 Sunday 16:30:13
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2025-04-24 Thursday 16:02:11
'''

import hashlib
import os
import requests
from bs4 import BeautifulSoup
import pandas as pd
from datetime import datetime
import argparse


def extract_href_from_tax_tree(tax_div):
    item = {}
    for a_tag in tax_div.find_all('a', href=True):
        rank, name = a_tag['href'].strip('/').split('/')
        item[rank.capitalize()] = name.capitalize()
    return item


def parse_genus_html(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    tax_trees = soup.find_all('div', class_='tax-tree')
    subtax = [extract_href_from_tax_tree(tax_tree) for tax_tree in tax_trees]
    return subtax


def fetch_web_content(url, cache_dir='lpsnCache'):
    # fetch web content from URL and save to cache if not already present
    os.makedirs(cache_dir, exist_ok=True)
    url_hash = hashlib.md5(url.encode('utf-8')).hexdigest()
    cache_file_path = os.path.join(cache_dir, f"{url_hash}.html")

    # check if the cache file exists
    if os.path.exists(cache_file_path):
        with open(cache_file_path, 'r') as file:
            return file.read()

    try:
        # fetch the web content from the URL
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        with open(cache_file_path, 'w') as file:
            file.write(response.text)
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Failure to fetch web content: {e}")
        return None


def get_taxtree_genus():
    taxSys = []
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for lab in list(letters):
        target_url = f"https://lpsn.dsmz.de/genus?page={lab}"
        html_content = fetch_web_content(target_url)
        data = parse_genus_html(html_content)
        taxSys.extend(data)
    lpsn = pd.DataFrame(taxSys)
    lpsn = lpsn.drop_duplicates()
    lpsn = lpsn.map(lambda x: 'Unclassified'
                    if isinstance(x, str) and '-no-' in x else x)
    return lpsn


def get_lpsn(outfile):
    if os.path.exists(outfile):
        lpsn = pd.read_csv(outfile)
    else:
        lpsn = get_taxtree_genus()
        lpsn.to_csv(outfile, index=False)
    return lpsn


mapheader = [('Domain(D)', 'Domain'),
             ('Phylum(P)', 'Phylum'),
             ('Class(C)', 'Class'),
             ('Order(O)', 'Order'),
             ('Family(F)', 'Family')]


def update_taxon(lpsn, infile, outfile, outmap=False):
    # update the taxon file by merging with the lpsn data
    taxon = pd.read_csv(infile)
    merged = pd.merge(taxon,
                      lpsn,
                      left_on='Genus(G)',
                      right_on='Genus',
                      how='left')
    if outmap:  # output the mapping file for check result
        merged.to_csv('map_taxon_lpsn.csv', index=False)

    # update the taxon file by lpsn data
    mask = merged['Genus'].notna()
    for t, l in mapheader:
        taxon.loc[mask, t] = merged.loc[mask, l]
    taxon.to_csv(outfile, index=False)

    # for row not find in LPSN
    nFail = (~mask).sum()
    if nFail > 0:
        print(f'there are {nFail} rows in {infile} not find in LPSN')
        taxon[~mask].to_csv(
            f'{os.path.splitext(infile)[0]}-lpsn-notfind.csv', index=False)


def parse_args():
    # the default LPSN file name is generated based on the current date
    current_date = datetime.now()
    lpsn_file = f'LPSN-GenusList-{current_date.strftime("%Y%m%d")}.csv'
    # set up the argument parser
    parser = argparse.ArgumentParser(
        description="Update taxon by LPSN with Genus Name")
    parser.add_argument(
        '-l',
        "--lpsn",
        type=str,
        default=lpsn_file,
        required=False,
        help="LPSN with Genus Name")
    parser.add_argument(
        '-o',
        "--outfile",
        type=str,
        default=None,
        required=False,
        help="Output new taxon file"
    )
    parser.add_argument(
        '-i',
        "--infile",
        type=str,
        required=False,
        help="Input taxon file, if not provided, output LPSN file only"
    )
    parser.add_argument(
        '-m',
        "--mapfile",
        action="store_true",
        help="Ouput mapfile map_taxon_lpsn.csv for check"
    )

    return parser.parse_args()


if __name__ == "__main__":
    # parse arguments
    args = parse_args()

    # get LPSN data
    lpsn = get_lpsn(args.lpsn)

    # update taxon file with LPSN data
    if args.infile:
        if not args.outfile:
            args.outfile = f'{os.path.splitext(args.infile)[0]}-lpsn.csv'
        update_taxon(lpsn, args.infile, args.outfile, args.mapfile)
