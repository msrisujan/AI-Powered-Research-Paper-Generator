from langchain_huggingface import HuggingFaceEndpoint
import os
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
from scholarly import scholarly
import arxiv
import re
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import time
import json
from fpdf import FPDF





rq = "How does artificial intelligence impact healthcare efficiency?"
keywords = "artificial intelligence, healthcare, efficiency"
ps = "Identify the specific impacts of AI on healthcare processes to improve efficiency."

sec_key = "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"
os.environ["HUGGINGFACEHUB_API_TOKEN"]=sec_key

repo_id="mistralai/Mistral-7B-Instruct-v0.3"
llm=HuggingFaceEndpoint(repo_id=repo_id,max_length=128,temperature=0.4,token=sec_key)

template="""
As an AI research assistant, your task is to refine and enhance the clarity of the following research inputs:

    Research Question: {research_question}
    Keywords: {keywords}
    Problem Statement: {problem_statement}

    Please analyze these inputs and provide a refined version that:
    1)can be passed as query to search for the existing papers related to the topic
    2)Provide just one sentence which is clear, concise, and improved.
    
just give refined research question as output with no explanation or headers or ask for assistance or note.
make sure no other text is included except the refined research question."""

refined_generation_prompt = PromptTemplate(
    input_variables=["research_question", "keywords", "problem_statement"],
    template=template)

refine_chain = refined_generation_prompt | llm

refined=refine_chain.invoke({
    "research_question": rq,
    "keywords": keywords,
    "problem_statement": ps
})

pattern = r"^\s*\n"
    # Use re.sub to replace empty lines with an empty string
refined = re.sub(pattern, "", refined, flags=re.MULTILINE)
print(refined)

results = []

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

#PubMed
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

#ACM Digital

def get_full_abstract(url):
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
                
                abstract = get_full_abstract(link)
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

# Example usage
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
}

search_query = refined[1:-1]
search_acm_library(search_query)

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

# Example usage
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
}

search_query = refined[1:-1]
search_springer_link(search_query)

repo_id="meta-llama/Meta-Llama-3-8B-Instruct"
summarizer=HuggingFaceEndpoint(repo_id=repo_id,max_length=128,temperature=0.4,token=sec_key)

template = """
        Based on the following information about a research paper, please provide:
        1. The main findings of the paper
        2. The key research questions addressed
        3. A brief summary of the paper

        Title: {title}
        Abstract: {abstract}

        Your response should be structured as follows:
        Findings: [Main findings of the paper]
        Research Questions: [Key research questions addressed]
        Summary: [A brief summary of the paper]
        
        no need to add other statements like "the main findings of paper are:", just seperate with commas.
        Make sure to follow the structure and do not generate any other text.
        """

summarize_prompt = PromptTemplate(
    input_variables=["title", "abstract"],
    template=template)

summary_chain = summarize_prompt | summarizer

def summarize_paper(paper):
    summary=summary_chain.invoke({
        "title": paper['Title'],
        "abstract": paper['Summary']
    })
    return summary

enhanced_results = []
for paper in results:
    summary = summarize_paper(paper)
    summary_dict = dict(item.split(": ", 1) for item in summary.split("\n") if ": " in item)
    paper.update(summary_dict)
    enhanced_results.append(paper)


def json_out(input_string):
    json_pattern = re.compile(r'{.*?}', re.DOTALL)
    json_match = json_pattern.search(input_string)
    json_string = json_match.group(0)
    json_data = json.loads(json_string)
    return json_data

template = """
    Based on the following research input and selected papers, generate a title and keywords for a research paper.

    Research Input:
    1)Refined Question: {refined_input}
    2)Keywords: {keywords}
    3)Problem Stamenet: {ps}
    

    Selected Papers:
    {paper_info}
    
        Your task is to:

    1. Create a Title:
        - The title should be clear, concise, and accurately reflect the main focus of the research
        - It should be engaging and capture the reader's interest
        - Aim for a length of 10-15 words
        - Consider including key variables or the main concept being studied
        - Avoid using abbreviations or jargon
        - If applicable, mention the study design (e.g., "A Randomized Controlled Trial of...")

    2. Generate Keywords:
        - Provide 5-7 keywords or short phrases
        - Keywords should reflect the main topics, concepts, or variables of the study
        - Include both general and specific terms related to the research
        - Consider including:
            * The main subject area
            * Key variables or concepts
            * Research methodology or study design
            * Population or sample characteristics (if relevant)
            * Geographical location (if relevant)

    Please structure your response as a JSON object with the following format:

    {{
        "Title": "Your generated title",
        "Keywords": "generated keywords",
    }}

    """

title_prompt = PromptTemplate(
    input_variables=["refined_input", "keywords", "ps", "paper_info"],
    template=template)
title_chain = title_prompt | summarizer
title_res = title_chain.invoke({
        "refined_input": refined,
        "keywords": keywords,
        "ps": ps,
        "paper_info": enhanced_results[:5]
       })
title_json = json_out(title_res)
title = title_json['Title']
keywords = title_json['Keywords']

template = """
        Generate an abstract for a research paper with the following title, keywords, Research Input and sample papers:

        Title: {title}
        Keywords: {keywords}

        Research Input:
        1)Refined Question: {refined_input}
        2)Problem Stamenet: {ps}
        
        Sample_papers:
        {paper_info}

        Your task is to create an abstract that includes the following components:

    1. Background: 
       - Briefly introduce the research topic and its importance
       - Provide context for the study

    2. Objective:
       - Clearly state the main research question or objective
       - If applicable, mention any hypotheses

    3. Methods:
       - Summarize the key aspects of the methodology
       - Include information on study design, participants, and main procedures

    4. Results:
       - Present the most important findings
       - Include key statistical results if applicable

    5. Conclusion:
       - Summarize the main conclusion(s) of the study
       - Briefly mention implications or applications of the findings

    Guidelines:
    - The abstract should be between 200 and 300 words
    - Use clear, concise language
    - Avoid jargon or unexplained acronyms
    - Do not include citations
    - Present information in a logical flow

    Please structure your response as a JSON object with the following format:

    {{
        "Abstract": "Your generated abstract"
    }}
        """

abstract_prompt = PromptTemplate(
    input_variables=["title", "keywords","refined_input", "ps", "paper_info"],
    template=template)

abstract_chain = abstract_prompt | summarizer

abstract_res = abstract_chain.invoke({
        "refined_input": refined,
        "keywords": keywords,
        "ps": ps,
        "paper_info": enhanced_results[:5],
        "title": title
       })
abstract_json = json_out(abstract_res)
abstract = abstract_json['Abstract']
print(abstract)

template = """
        Generate an introduction for a research paper with the following title, abstract, research input and sample papers:

        Title: {title}
        Abstract: {abstract}

        Research Input:
        1)Refined Question: {refined_input}
        2)Problem Stamenet: {ps}
        
        Sample_papers:
        {paper_info}

         Your task is to create an introduction that includes the following components:

    1. Background:
       - Provide context for the research topic
       - Explain the importance and relevance of the study
       - Briefly introduce key concepts or theories related to the research

    2. Problem Statement:
       - Clearly articulate the problem or gap in knowledge that the research addresses
       - Explain why this problem is significant and worthy of study

    3. Literature Review:
       - Summarize relevant previous research
       - Identify gaps or contradictions in existing literature
       - Explain how your study relates to or builds upon previous work

    4. Research Question and Objectives:
       - Clearly state the main research question(s) or objective(s)
       - If applicable, present any hypotheses

    5. Significance of the Study:
       - Explain the potential contributions of your research
       - Discuss possible implications for theory, practice, or policy

    6. Overview of Methodology:
       - Briefly introduce the research approach or design
       - Do not go into detailed methods (save that for the methodology section)

    7. Structure of the Paper:
       - Provide a brief overview of how the rest of the paper is organized

    Guidelines:
    - The introduction should be between 500 and 750 words
    - Use clear, academic language
    - Cite relevant literature (use placeholder citations like [Author, Year])
    - Ensure a logical flow of ideas
    - End with a clear transition to the next section of the paper

    Please structure your response as a JSON object with the following format:

    {{
        "Introduction": "Your generated introduction"
    }}

    """

intro_prompt = PromptTemplate(
    input_variables=["title", "abstract","refined_input", "ps", "paper_info"],
    template=template)

intro_chain = intro_prompt | summarizer

intro_res = intro_chain.invoke({
        "refined_input": refined,
        "abstract": abstract,
        "ps": ps,
        "paper_info": enhanced_results[:5],
        "title": title
       })
# print(intro_res)
intro_json = json_out(intro_res)
introduction = intro_json['Introduction']
print(introduction)

template = """
    Create a literature review based on the following:

    Research Question: {refined_input}
    Keywords: {keywords}
    
    Selected Papers:
    {paper_info}

    Include:
    1. Brief introduction
    2. Key themes and findings
    3. Gaps in current research
    4. How your study fits in
    5. Short summary

    Guidelines:
    - Aim for 800-1000 words.
    - Cite all referenced studies in order of selected papers as numbers
    - Be objective in presenting previous research.
    - Ensure logical flow between topics.

    Please structure your response as a JSON object with the following format:

    {{
        "Literature Review": "Your generated literature review",
    }}
    Ensure that there is only one key ('Literature Review') and one value (which should not be json) in which the value contains all the paragraphs of the Literature Review section.

    """

literature_prompt = PromptTemplate(
    input_variables=["refined_input", "keywords", "paper_info"],
    template=template)

literature_chain = literature_prompt | summarizer

literature_res = literature_chain.invoke({
        "refined_input": refined,
        "keywords": keywords,
        "paper_info": enhanced_results[:5]
       })
# print(literature_res)
literature_json = json_out(literature_res)
literature = literature_json['Literature Review']
print(literature)

template = """
Generate a methodology section for a research paper based on the following research input and other sections:

Research Input:
1)Refined Question: {refined_input}
2)Problem Stamenet: {ps}


Title: {title}
Abstract: {abstract}
Introduction: {intro}
literature: {literature}

Structure the methodology as follows:

1. Research Design:
    - Briefly describe and justify the chosen research design.

2. Participants:
    - Describe the sample, including size, demographics, and selection criteria.

3. Data Collection:
    - Explain the methods and tools used to collect data.
    - Describe any instruments or measures in detail.

4. Procedure:
    - Outline the step-by-step process of conducting the study.

5. Data Analysis:
    - Describe the analytical methods used, including any statistical tests.

6. Ethical Considerations:
    - Briefly mention ethical approvals and considerations.

Guidelines:
- Aim for 700-800 words.
- Use clear, precise language.
- Provide enough detail for replication.
- Justify methodological choices where necessary.

Please structure your response as a JSON object with the following format:

{{
    "Methodology": "Your generated methodology",
}}

Ensure that there is only one key ('Methodology') and one value (which should not be json) in which the value contains all the paragraphs of the methodology section.

"""

methodology_prompt = PromptTemplate(
    input_variables=["refined_input", "ps", "title","abstract", "intro", "literature"],
    template=template)

methodology_chain = methodology_prompt | summarizer

methodology_res = methodology_chain.invoke({
        "refined_input": refined,
        "ps": ps,
        "title": title,
        "abstract": abstract,
        "intro": introduction,
        "literature": literature
       })
# print(methodology_res)
methodology_json = json_out(methodology_res)
# print(methodology_json)
methodology = methodology_json['Methodology']
print(methodology)

template = """
    Generate a detailed results section for a research paper based on the following information:

    Research Question: {refined_input}
    abstract: {abstract}
    Methodology Summary: {methodology}

    Your task is to create a comprehensive results section that presents the findings of the study. Remember, this section should focus on reporting the results without interpretation. Structure your response to include the following components:

    1. Overview of Data Analysis:
       - Briefly restate the main data analysis procedures used
       - Mention any data preprocessing or cleaning steps

    2. Descriptive Statistics:
       - Present key descriptive statistics relevant to the research question
       - Include measures of central tendency and variability where appropriate
       - Consider using bullet points or a table format for clarity

    3. Main Findings:
       - Present the primary results that directly address the research question
       - Organize findings logically, possibly in order of importance or as they relate to research sub-questions
       - Include specific data points, statistical test results, or qualitative findings as appropriate

    4. Secondary or Exploratory Findings:
       - Present any additional results that may be relevant but were not the primary focus
       - Explain how these findings relate to the main research question

    5. Visual Representations:
       - Suggest appropriate figures, charts, or tables to visually represent key findings
       - Provide brief descriptions of what each visual representation shows

    6. Statistical Significance:
       - Report p-values, confidence intervals, or other measures of statistical significance where relevant
       - Clearly state whether results support or refute the research hypotheses (if applicable)

    7. Unexpected Results:
       - Mention any unexpected or surprising findings
       - Avoid explaining or interpreting these results (save that for the discussion section)

    Use neutral, objective language to present the findings without bias or interpretation. If specific numerical results are not available, provide placeholder values or ranges that would be typical for this kind of study.
    
     Please structure your response in json format:
    {{
    "Results": "Your generated results here. Include all paragraphs in this single string value."
    }}
    Ensure that there is only one key ('Results') and one value (which should not be json) in which the value contains all the paragraphs of the results section.
    """

results_prompt = PromptTemplate(
    input_variables=["refined_input", "abstract", "methodology"],
    template=template)

results_chain = results_prompt | summarizer

results_res = results_chain.invoke({
        "refined_input": refined,
        "abstract": abstract,
        "methodology": methodology
       })
results_json = json_out(results_res)
# print(results_json)
results = results_json['Results']
print(results)

template = """
    Based on the following information about a research study, generate a comprehensive conclusion section.

    Title: {title}
    Abstract: {abstract}
    Research Question: {refined_input}
    Methodology: {methodology}
    Results: {results}

    Your task is to create a conclusion that includes the following components:

    1. Restatement of Research Purpose:
       - Briefly restate the main research question or objective
       - Remind the reader of the importance of the study

    2. Summary of Key Findings:
       - Concisely summarize the main results of the study
       - Relate these findings directly to the research question
       - Highlight any unexpected or particularly significant results

    3. Interpretation of Results:
       - Explain what the results mean in the context of your field
       - Discuss how the findings support or challenge existing theories or previous research

    4. Implications:
       - Discuss the theoretical implications of your findings
       - Explain practical applications or real-world relevance of your results
       - If applicable, suggest policy implications

    5. Limitations:
       - Acknowledge any limitations or weaknesses of the study
       - Explain how these limitations might affect the interpretation of the results

    6. Future Research Directions:
       - Suggest areas for future research based on your findings
       - Identify questions that remain unanswered or new questions raised by your study

    7. Concluding Statement:
       - Provide a strong final statement that encapsulates the main message of your study
       - Emphasize the contribution of your research to the field

    Guidelines:
    - The conclusion should be between 400 and 600 words
    - Use clear, academic language
    - Avoid introducing new information not previously discussed in the paper
    - Ensure a logical flow of ideas
    - End on a strong, memorable note

    Please structure your response as a JSON object with the following format:

    {{
        "Conclusion": "Your generated conclusion"
    }}
    """

conclusion_prompt = PromptTemplate(
    input_variables=["refined_input", "abstract", "methodology"],
    template=template)

conclusion_chain = conclusion_prompt | summarizer

conclusion_res = conclusion_chain.invoke({
        "refined_input": refined,
        "abstract": abstract,
        "methodology": methodology,
        "title": title,
        "results": results
       })
conclusion_json = json_out(conclusion_res)
# print(conclusion_json)
conclusion = conclusion_json['Conclusion']
print(conclusion)


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
        self.set_font('Arial', '', 12)
        for ref in references:
            title = ref['Title']
            link = ref['Link']
            self.set_text_color(0, 0, 255)  # Blue color for links
            self.cell(0, 10, title, ln=True, link=link)
            self.ln()
        self.set_text_color(0, 0, 0)


# Create PDF
pdf = PDF()

pdf.add_page()
pdf.chapter_title(title)
pdf.chapter_body(abstract)
pdf.chapter_title("Introduction")
pdf.chapter_body(introduction)
pdf.chapter_title("Literature Review")
pdf.chapter_body(literature)
pdf.chapter_title("Methodology")
pdf.chapter_body(methodology)
pdf.chapter_title("Results")
pdf.chapter_body(results)
pdf.chapter_title("Conclusion")
pdf.chapter_body(conclusion)
pdf.add_reference(enhanced_results)

# Save PDF
pdf_output_path = 'research_paper.pdf'
pdf.output(pdf_output_path)

if os.path.exists(pdf_output_path):
    print(f"PDF saved successfully at {pdf_output_path}")
else:
    print(f"Failed to save PDF at {pdf_output_path}")