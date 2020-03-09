########################################################################################################################
# 
# Lookup CAS number in both EPA DSSTox file and via CompTox dashboard.
# 
########################################################################################################################

import logging
import requests
import re
import pandas as pd

# logging.getLogger('requests').setLevel('WARNING')

from lxml import html
from rdkit import Chem

# from project_utils import *

########################################################################################################################

# Config...

dashboard_url_prefix = "https://comptox.epa.gov"
dashboard_url = dashboard_url_prefix + "/dashboard/dsstoxdb/results?utf8=✓&search={cas}"

########################################################################################################################

logger = logging.getLogger(__name__) # Elisabet Gregori: substitute import make_logger by logging 

########################################################################################################################
 
# HTML to link to CompTox Dashboard for a CAS number...

comptox_link = lambda x: '<a target="_blank" href="https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=✓&search={cas}">link</a>'.format(cas=x)

##########

def get_conntab(doc):

    """
    Extract URL from lxml document version of Dashboard page and retrieve molfile.
    """

    molfile_url = dashboard_url_prefix + doc.xpath('//a[starts-with(@href, "/dashboard/dsstoxdb/download_mol_file")]/@href')[0]

    reponse = requests.get(molfile_url)

    if reponse.status_code == 200:

        conntab = reponse.content.decode('utf-8')

    else:

        conntab = None

    return conntab

##########

def get_related(doc):

    """
    Get 'related chemicals' from lxml document version of Dashboard page.
    """

    try:

        li = doc.xpath('//div[@id="related-chemcials-container"]/div/li')[0]

        cas_new = li.xpath('./div[@class="chem-casrn"]/p/text()')[0]

        a = li.xpath('./div[@class="chem-name"]/p/a')[0]

        name_new, url_new = a.text.strip(), a.attrib['href']

    except (IndexError, ValueError) as error:

        name_new, url_new = None, None

    return name_new, url_new

##########

def comptox_lookup(cas, flag_related=False):

    """
    Get CompTox Dashboard page for a CAS number and extract structural information, identifiers etc.
    """

    columns = ['cas', 'name', 'inchi', 'inchikey', 'smiles', 'conntab', 'synonyms']

    url = dashboard_url.format(cas=cas)

    ######

    # Retrieve document...

    while True:

        response = requests.get(url)

        assert response.status_code == 200

        doc = html.fromstring(response.content.decode('utf-8'))

        message = ''.join(doc.xpath('''
                                  //div[@id="chemical-searched-by"]/p/text()
                                | //p[@id="advancedSearchedFor"]/text()
                                | //h1[text()="Search Results"]/following-sibling::p/text()
                                | //h1[text()="Search Results"]/following-sibling::div[2]/div/p/text()
                            ''')).strip()

        if re.match("^Searched by Deleted CAS-RN", message):

            cas_found, dtxsid = [x.strip() for x in doc.xpath('//h1[@id="chem-label"]/small[@id="casrn-subtitle"]/text()')[0].split('|')]

            logger.warning("CAS number '{cas}': deleted; CAS number actually found = '{cas_found}'.".format(cas=cas, cas_found=cas_found))

        elif re.match("^Searched by Synonym: Found [1-9][0-9]* results", message):

            hrefs = list(set(doc.xpath('//div[@id="multiple-results-view"]//div[contains(@class, "chem-dtxsid")]/p/a/@href')))

            if len(hrefs) == 1:

                url = dashboard_url_prefix + hrefs[0]

                logger.warning("CAS number '{cas}': multiple results but only one link, so will attempt to follow ('{url}').".format(cas=cas, url=url))

                continue

            else:

                logger.warning("CAS number '{cas}': multiple results ({n_hrefs}).".format(cas=cas, n_hrefs=len(hrefs)))

                return None

        elif re.match("^(?:Searched by Synonym: Found 0 results|No results found for ''.)", message):

            logger.debug("CAS number '{cas}': no hits.".format(cas=cas))

            return None

        break

    ######

    # Extract name...

    try:

        name = doc.xpath('//h1[@id="chem-label"]/@data-label')[0]

    except (IndexError) as error:

        logger.warning("CAS number '{cas}': problem extracting name.".format(cas=cas))

        name = None

    ######

    # Extract synonyms...

    synonyms = [y[0].replace('\\n', '') for y in (x.xpath('./text() | ./*/text()') for x in doc.xpath('//div[@id="synonyms"]/div/table/tr/td')) if y]

    ######

    # Extract structures...

    try:

        divs = doc.xpath('//div[@id="structural-identifiers-panel"]/div[2]/div/div')

        iupac_name, smiles, inchi, inchikey = [''.join(div.xpath('./div/p/text()')).strip().split('\n')[0] for div in divs[:4]]

        conntab = get_conntab(doc)

        if conntab and (not smiles or smiles == 'Not Found'):

            logger.warning("CAS number '{cas}' (name = '{name}'): will regenerate SMILES from conntab.".format(cas=cas, name=name))

            mol = Chem.MolFromMolBlock(conntab)

            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    except (IndexError, ValueError) as error: # No structure available.

        logger.warning("CAS number '{cas}' (name = '{name}'): no structure available.".format(cas=cas, name=name))

        if flag_related:

            name_new, url_new = get_related(doc)

            logger.warning("\trelated chemical suggested: cas = '{cas_new}', name = '{name_new}', url = '{url_new}'.".format(cas=cas, name=name, cas_new=cas_new, name_new=name_new, url_new=url_new))

        return None # inchi, inchikey, smiles, conntab = None, None, None, None

    ######

    # Finished...

    return pd.Series([cas, name, inchi, inchikey, smiles, conntab, '<br>'.join(synonyms)], index=columns)

# comptox_lookup()

########################################################################################################################
# End
########################################################################################################################
