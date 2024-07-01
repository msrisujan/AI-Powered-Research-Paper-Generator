import arxiv
import re
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import time
import json
from fpdf import FPDF



headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
}

results = []

def collect(refined):

    # arxiv
    client = arxiv.Client()

    search = arxiv.Search(
    query = refined[1:-1],
    max_results = 5,
    sort_by = arxiv.SortCriterion.Relevance,
    sort_order = arxiv.SortOrder.Descending
    )

    for result in client.results(search):
        results.append({
            'Title': result.title,
            'Source': 'arXiv',
            'Link': result.entry_id,
            'Findings': 'N/A',  
            'Research Questions': 'N/A', 
            'Summary': result.summary
        })

    Entrez.email = "srisujanbannu@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=refined[1:-1], retmax=10)
    record = Entrez.read(handle)
    id_list = record["IdList"]
    id_list = id_list[:5]
    papers = []
    for pmid in id_list:
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        fetch_record = fetch_handle.read()
        title = ''
        abstract = ''
        capture_title = False
        capture_abstract = False

        lines = fetch_record.split('\n')
        for line in lines:
            if line.startswith("TI  - "):
                capture_title = True
                title = line[6:]
            elif capture_title and line.startswith("      "):  
                title += ' ' + line.strip()
            else:
                capture_title = False

            if line.startswith("AB  - "):
                capture_abstract = True
                abstract = line[6:]
            elif capture_abstract and line.startswith("      "):  
                abstract += ' ' + line.strip()
            else:
                capture_abstract = False

        results.append({
            'Title': title,
            'Source': 'PubMed',
            'Link': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            'Findings': 'N/A',  
            'Research Questions': 'N/A',  
            'Summary': abstract
        })

    search_acm_library(refined[1:-1])
    search_springer_link(refined[1:-1])

    return results


def get_full_abstract_acm(url):
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        abstract_elem = soup.find('div', {'role': 'paragraph'})
        if abstract_elem:
            return abstract_elem.text.strip()
    return "Abstract not available"

def search_acm_library(query, max_results=5):
    base_url = "https://dl.acm.org/action/doSearch"
    params = {
        "AllField": query,
        "startPage": "0",
        "pageSize": str(max_results)
    }
    
    response = requests.get(base_url, params=params, headers=headers)
    
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        res = soup.find_all('div', class_='issue-item__content')
        
        for i, result in enumerate(res[:max_results], 1):
            title_elem = result.find('h5', class_='issue-item__title')
            
            if title_elem:
                title = title_elem.text.strip()
                link = "https://dl.acm.org" + title_elem.find('a')['href']
                
                abstract = get_full_abstract_acm(link)
                results.append({
                    'Title': title,
                    'Source': 'ACM Digital Library',
                    'Link': link,
                    'Findings': 'N/A',  
                    'Research Questions': 'N/A',  
                    'Summary': abstract
                })
                
                # Add a small delay to avoid overwhelming the server
                time.sleep(1)
            else:
                print(f"Paper {i}: Unable to extract title")
                print("-" * 50)
    else:
        print(f"Error: Unable to fetch results. Status code: {response.status_code}")


#springer link
def get_full_abstract(url):
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        abstract_elem = soup.find('div', class_='c-article-section__content', id='Abs1-content')
        if abstract_elem:
            return abstract_elem.text.strip()
    return "Abstract not available"

def search_springer_link(query, max_results=5):
    base_url = "https://link.springer.com/search"
    params = {
        "query": query,
        "facet-content-type": "Article"
    }
    
    response = requests.get(base_url, params=params, headers=headers)
    
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        res = soup.find_all('li', class_='app-card-open')
        
        for i, result in enumerate(res[:max_results], 1):
            title_elem = result.find('h3')
            
            if title_elem:
                title = title_elem.text.strip()
                link = "https://link.springer.com" + title_elem.find('a')['href']
                
                abstract = get_full_abstract(link)
                results.append({
                    'Title': title,
                    'Source': 'Springer Link',
                    'Link': link,
                    'Findings': 'N/A',  
                    'Research Questions': 'N/A',  
                    'Summary': abstract
                })
                
                # Add a small delay to avoid overwhelming the server
                time.sleep(1)
            else:
                print(f"Paper {i}: Unable to extract title")
                print("-" * 50)
    else:
        print(f"Error: Unable to fetch results. Status code: {response.status_code}")


class PDF(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, 'Research Paper', 0, 1, 'C')

    def chapter_title(self, title):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(5)

    def chapter_body(self, body):
        self.set_font('Arial', '', 12)
        self.multi_cell(0, 10, body)
        self.ln()

    def add_reference(self, references):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, 'References', 0, 1, 'L')
        self.ln(5)
        self.set_font('Arial', '', 10)
        for ref in references:
            title = ref['Title']
            link = ref.get('Link', '')  
            
            # Add the title
            self.set_text_color(0, 0, 0)  
            self.multi_cell(0, 5, title)
            
            # Add the link if available
            if link:
                self.set_text_color(0, 0, 255) 
                self.multi_cell(0, 5, link)
                self.link(self.get_x(), self.get_y() - 5, self.get_string_width(link), 5, link)
            
            self.ln(5)  
        
        self.set_text_color(0, 0, 0) 
