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
import streamlit as st
from dotenv import load_dotenv
import pandas as pd

from prompts import *
from collect_papers import collect
from collect_papers import PDF

load_dotenv()

sec_key = os.getenv("HUGGINGFACEHUB_API_TOKEN")

os.environ["HUGGINGFACEHUB_API_TOKEN"] = sec_key

repo_id = "mistralai/Mistral-7B-Instruct-v0.3"
llm = HuggingFaceEndpoint(repo_id=repo_id, max_length=128, temperature=0.4, token=sec_key)

st.title("AI Research Assistant")

# Taking inputs from the user
rq = st.text_input("Research Question", "How does artificial intelligence impact healthcare efficiency?")
keywords = st.text_input("Keywords", "artificial intelligence, healthcare, efficiency")
ps = st.text_input("Problem Statement", "Identify the specific impacts of AI on healthcare processes to improve efficiency.")

refined_generation_prompt = PromptTemplate(
    input_variables=["research_question", "keywords", "problem_statement"],
    template=refine_template
)

refine_chain = refined_generation_prompt | llm

def refine_research_question(rq, keywords, ps):
    refined = refine_chain.invoke({
        "research_question": rq,
        "keywords": keywords,
        "problem_statement": ps
    })
    pattern = r"^\s*\n"
    refined = re.sub(pattern, "", refined, flags=re.MULTILINE)
    return refined

# State to store refined question
if 'refined_rq' not in st.session_state:
    st.session_state.refined_rq = None

if 'continue_process' not in st.session_state:
    st.session_state.continue_process = False

# Button to refine research question
if st.button('Generate Refined Research Question'):
    with st.spinner('Refining research question...'):
        st.session_state.refined_rq = refine_research_question(rq, keywords, ps)

# If refined research question is available, show buttons to continue or regenerate
if st.session_state.refined_rq and not st.session_state.continue_process:
    st.write("Refined Research Question:", st.session_state.refined_rq)
    col1, col2 = st.columns(2)
    with col1:
        if st.button('Continue'):
            st.session_state.continue_process = True
            st.experimental_rerun()
    with col2:
        if st.button('Regenerate'):
            with st.spinner('Refining research question...'):
                st.session_state.refined_rq = refine_research_question(st.session_state.refined_rq, keywords, ps)
                st.experimental_rerun()

repo_id = "meta-llama/Meta-Llama-3-8B-Instruct"
summarizer = HuggingFaceEndpoint(repo_id=repo_id, max_length=128, temperature=0.4, token=sec_key)

def summarize_paper(paper):
    summary = summary_chain.invoke({
        "title": paper['Title'],
        "abstract": paper['Summary']
    })
    return summary

summarize_prompt = PromptTemplate(
            input_variables=["title", "abstract"],
            template=enhance_template
        )
summary_chain = summarize_prompt | summarizer

if 'enhanced_results' not in st.session_state:
    st.session_state.enhanced_results = None

def collect_and_enhance_papers():
    with st.spinner('Collecting papers from each source...'):
        results = collect(st.session_state.refined_rq)

    st.write("Papers collected...")

    with st.spinner('Enhancing the collected results...'):
        
        enhanced_results = []
        for paper in results:
            summary = summarize_paper(paper)
            summary_dict = dict(item.split(": ", 1) for item in summary.split("\n") if ": " in item)
            paper.update(summary_dict)
            enhanced_results.append(paper)
    
    return enhanced_results

# Add this near the top of your script, with other session state initializations
if 'show_enhanced_results' not in st.session_state:
    st.session_state.show_enhanced_results = True
if 'generate_paper' not in st.session_state:
    st.session_state.generate_paper = False
if 'selected_papers' not in st.session_state:
    st.session_state.selected_papers = []

# Replace the existing code from "if 'continue_process' in st.session_state..." onwards with this:
def json_out(input_string):
    json_pattern = re.compile(r'{.*?}', re.DOTALL)
    json_match = json_pattern.search(input_string)
    json_string = json_match.group(0)
    json_data = json.loads(json_string)
    return json_data

if 'continue_process' in st.session_state and st.session_state.continue_process:
    if st.session_state.enhanced_results is None:
        st.session_state.enhanced_results = collect_and_enhance_papers()

    enhanced_results = st.session_state.enhanced_results

    if st.session_state.show_enhanced_results:
        st.write('## Select up to 5 papers')
        selected_papers = []
        if enhanced_results:
            for idx, paper in enumerate(enhanced_results):
                st.write(f"### Paper {idx + 1}")
                st.write("Title:", paper['Title'])
                st.write("Source:", paper['Source'])
                st.write("Findings:", paper['Findings'])
                st.write("Research Questions:", paper['Research Questions'])
                st.write("Summary:", paper['Summary'])
                selected = st.checkbox("Select", key=f"paper_{idx}")
                if selected:
                    selected_papers.append(paper)
                if len(selected_papers) > 5:
                    st.warning("You can select up to 5 papers only.")
                    selected_papers.pop()

        if st.button('Select these Papers'):
            st.session_state.show_enhanced_results = False
            st.session_state.selected_papers = selected_papers
            st.experimental_rerun()

    else:  # Show table of selected papers
        st.write("### Table of Selected Papers")
        papers_df = pd.DataFrame(st.session_state.selected_papers)
        st.dataframe(papers_df[['Title', 'Source', 'Findings', 'Research Questions', 'Summary']])
        col1, col2 = st.columns(2)
        with col1:
            if st.button('Generate Paper'):
                st.session_state.generate_paper = True
        with col2:
            if st.button('Reselect Papers'):
                st.session_state.show_enhanced_results = True
                st.experimental_rerun()

    if st.session_state.generate_paper:
        with st.spinner('Generating Title and Keywords..'):
            title_prompt = PromptTemplate(
                input_variables=["refined_input", "keywords", "ps", "paper_info"],
                template=title_template)
            title_chain = title_prompt | summarizer
            title_res = title_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "keywords": keywords,
                    "ps": ps,
                    "paper_info": st.session_state.selected_papers
                })
            title_json = json_out(title_res)
            title = title_json['Title']
            keywords = title_json['Keywords']
        st.write('Title and Keywords Generated..')
        with st.spinner('Generating Abstract...'):
            abstract_prompt = PromptTemplate(
                input_variables=["title", "keywords","refined_input", "ps", "paper_info"],
                template=abstract_template)

            abstract_chain = abstract_prompt | summarizer

            abstract_res = abstract_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "keywords": keywords,
                    "ps": ps,
                    "paper_info": st.session_state.selected_papers,
                    "title": title
                })
            try:
                abstract_json = json_out(abstract_res)
                abstract = abstract_json['Abstract']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                abstract = abstract_res
        st.write('Abstract Generated..')
        with st.spinner('Generating Introduction...'):
            intro_prompt = PromptTemplate(
                input_variables=["title", "abstract","refined_input", "ps", "paper_info"],
                template=intro_template)

            intro_chain = intro_prompt | summarizer

            intro_res = intro_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "abstract": abstract,
                    "ps": ps,
                    "paper_info": st.session_state.selected_papers,
                    "title": title
                })
            # print(intro_res)
            try:
                intro_json = json_out(intro_res)
                introduction = intro_json['Introduction']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                introduction = intro_res
        st.write('Introduction generated...')
        with st.spinner('Generating Literature Review...'):
            literature_prompt = PromptTemplate(
                input_variables=["refined_input", "keywords", "paper_info"],
                template=literature_template)

            literature_chain = literature_prompt | summarizer

            literature_res = literature_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "keywords": keywords,
                    "paper_info": st.session_state.selected_papers
                })
            # print(literature_res)
            try:
                literature_json = json_out(literature_res)
                literature = literature_json['Literature Review']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                literature = literature_res
        st.write('Literature Review Generated...')
        with st.spinner('Generating Methodology...'):
            methodology_prompt = PromptTemplate(
                input_variables=["refined_input", "ps", "title","abstract", "intro", "literature"],
                template=methodology_template)

            methodology_chain = methodology_prompt | summarizer

            methodology_res = methodology_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "ps": ps,
                    "title": title,
                    "abstract": abstract,
                    "intro": introduction,
                    "literature": literature
                })
            # print(methodology_res)
            try:
                methodology_json = json_out(methodology_res)
                # print(methodology_json)
                methodology = methodology_json['Methodology']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                methodology = methodology_res
        st.write('Methodology Generated...')
        with st.spinner('Generating Results...'):
            results_prompt = PromptTemplate(
                input_variables=["refined_input", "abstract", "methodology"],
                template=results_template)

            results_chain = results_prompt | summarizer

            results_res = results_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "abstract": abstract,
                    "methodology": methodology
                })
            try:
                results_json = json_out(results_res)
                # print(results_json)
                results = results_json['Results']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                results = results_res
        st.write('Results Generated...')
        with st.spinner('Generating Conclusion...'):
            conclusion_prompt = PromptTemplate(
                input_variables=["refined_input", "abstract", "methodology"],
                template=conclusion_template)

            conclusion_chain = conclusion_prompt | summarizer

            conclusion_res = conclusion_chain.invoke({
                    "refined_input": st.session_state.refined_rq,
                    "abstract": abstract,
                    "methodology": methodology,
                    "title": title,
                    "results": results
                })
            try:
                conclusion_json = json_out(conclusion_res)
                # print(conclusion_json)
                conclusion = conclusion_json['Conclusion']
            except Exception as e:
                st.warning(f"Error parsing JSON: {str(e)}. Using raw output.")
                conclusion = conclusion_res
        st.write('Conclusion Generated')
        pdf = PDF()

        pdf.add_page()
        pdf.set_font("Arial", "B", 16)
        pdf.multi_cell(0, 10, title, 0, 'C')
        pdf.set_font("Arial", "I", 12)
        pdf.multi_cell(0, 10, f"Keywords: {keywords}", 0, 'C')
        pdf.ln(10)
        pdf.chapter_title("Abstract")
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
        pdf.add_reference(st.session_state.selected_papers)

        # Save PDF
        pdf_output_path = 'research_paper.pdf'
        pdf.output(pdf_output_path)

        with open(pdf_output_path, "rb") as pdf_file:
            pdf_bytes = pdf_file.read()

        # Display the download button
        st.download_button(
            label="Download Research Paper",
            data=pdf_bytes,
            file_name="research_paper.pdf",
            mime="application/pdf"
        )

        # Display a success message
        st.success("Research paper generated successfully! Click the button above to download.")

        # Reset the generate_paper state
        st.session_state.generate_paper = False


        