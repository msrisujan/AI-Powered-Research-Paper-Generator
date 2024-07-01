

refine_template="""
As an AI research assistant, your task is to refine and enhance the clarity of the following research inputs:

    Research Question: {research_question}
    Keywords: {keywords}
    Problem Statement: {problem_statement}

    Please analyze these inputs and provide a refined version that:
    1)can be passed as query to search for the existing papers related to the topic
    2)Provide just one sentence which is clear, concise, and improved.
    
just give refined research question as output with no explanation or headers or ask for assistance or note.
make sure no other text is included except the refined research question."""

enhance_template = """
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

title_template = """
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

abstract_template = """
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

intro_template = """
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

literature_template = """
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

methodology_template = """
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

results_template = """
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

conclusion_template = """
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