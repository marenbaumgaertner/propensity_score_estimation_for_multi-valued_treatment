import random
import pandas as pd
import numpy as np
import re
import warnings
from lxml import html, etree

import requests
from bs4 import BeautifulSoup
import time
import os
from selenium import webdriver
from selenium.webdriver.common.keys import Keys

#from selenium.webdriver.chromium import webdriver
#from selenium.webdriver.edge.service import Service
#from selenium.webdriver import EdgeOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

options = webdriver.ChromeOptions()
options.add_argument('user-agent=')
options.add_argument('--disable-blink-features=AutomationControlled')

# define URL
url = "https://topepo.github.io/caret/train-models-by-tag.html#Supports_Class_Probabilities"

# Create a Selenium WebDriver using Chrome
#driver = webdriver.Chrome(executable_path='path_to_chromedriver')
#driver = webdriver.Chrome("C:/Users/m-bau/miniconda3/Lib/site-packages/selenium/webdriver/chrome/chromedriver.exe")
#driver.get(url)
#
#subheader_elements = driver.find_element(By.XPATH, '//*[@id="supports-class-probabilities"]/p[2]/strong').text
#subheaders = [element.text for element in subheader_elements]

# Send an HTTP GET request to the website
response = requests.get(url)


# Check if the request was successful
if response.status_code == 200:
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')

    subheaders = []
    hyperparams = []
    j = 7

    for i in range(2, 698):  # Loop until 698
        subheader_element = soup.select(f'#supports-class-probabilities p:nth-of-type({i}) strong')
        if subheader_element:
            subheaders.append((i, subheader_element[0].get_text()))
            if soup.select(f'#supports-class-probabilities p:nth-of-type({i+2})')[0].get_text() == 'No tuning parameters for this model':
                hyperparams.append((i, 'no hyperparameters'))
                continue
            while not soup.select(f'#supports-class-probabilities > ul:nth-child({j})'):
                #print(i, j)
                j += 1
            k=1
            while soup.select(f'#supports-class-probabilities > ul:nth-child({j}) > li:nth-child({k})'):
                hyperparams_element = soup.select(f'#supports-class-probabilities > ul:nth-child({j}) > li:nth-child({k})')
                hyperparams.append((i, hyperparams_element[0].get_text()))
                #print(i, j, k)
                k += 1

            j += 1
            print(i,j, subheader_element[0].get_text())


else:
    print("Failed to retrieve the web page. Status code:", response.status_code)

# Create a DataFrame from the subheaders
df = pd.DataFrame({'Subheaders': subheaders})

# Export the DataFrame to a CSV file
df.to_csv('subheaders.csv', index=False)

# Create dictionaries to store the data
data = {"Subheader": [], "Hyperparams": []}

# Extract subheaders and corresponding hyperparameters
for i, subheader in subheaders:
    data["Subheader"].append(subheader)
    hyperparam_list = [hyperparam for j, hyperparam in hyperparams if i == j]
    data["Hyperparams"].append(", ".join(hyperparam_list))

# Create a data frame and set 'i' as the index
df = pd.DataFrame(data)
#df.set_index(pd.Index([i for i, _ in subheaders]), inplace=True)
# Export the DataFrame to a CSV file
df.to_csv('subheaders_new.csv', index=False)

# Display the data frame
print(df)


# Check if the request was successful
if response.status_code == 200:
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')

    subheaders_binary_only = []

    for i in range(0, 137):  # Loop until 698 with a step of 4
        subheader_element = soup.select(f'#two-class-only p:nth-of-type({i}) strong')
        if subheader_element:
            subheaders_binary_only.append(subheader_element[0].get_text())

else:
    print("Failed to retrieve the web page. Status code:", response.status_code)

df['binary_only'] = df['Subheader'].apply(lambda x: 1 if x in subheaders_binary_only else 0)
df.to_csv('subheaders.csv', index=False)


# Check if the request was successful
if response.status_code == 200:
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')

    subheaders_tree = []

    for i in range(0, 139):  # Loop until 698 with a step of 4
        subheader_element = soup.select(f'#tree-based-model p:nth-of-type({i}) strong')
        if subheader_element:
            subheaders_tree.append(subheader_element[0].get_text())

else:
    print("Failed to retrieve the web page. Status code:", response.status_code)

df['tree_based'] = df['Subheader'].apply(lambda x: 1 if x in subheaders_tree else 0)
df.to_csv('subheaders.csv', index=False)


len(df['Subheaders'].unique())


########################################################################################################################

//*[@id="supports-class-probabilities"]/p[3]
//*[@id="supports-class-probabilities"]/p[4]
//*[@id="supports-class-probabilities"]/ul[1]/li[1]/code/text()
//*[@id="supports-class-probabilities"]/ul[2]/li[1]/code/text()

for ul_html in hyperparams_elements:
    # Parse the <ul> element using BeautifulSoup
    hyperparams_soup = BeautifulSoup(hyperparams_element, 'html.parser')

    # Find all <li> elements within the <ul>
    li_elements = hyperparams_soup.find_all('li')

    for li_element in li_elements:
        # Extract and print the text
        text = li_element.get_text()
        print(text)

hyperparams_element = soup.find(f'// *[ @ id = "supports-class-probabilities"] / ul[2]')
element = soup.find('div', class_='supports-class-probabilities').find_all('li')[0]

# Iterate over <li> elements within the <ul>
for li_element in hyperparams_elements.find_all('li'):
    code_element = li_element.find('code')
    text = code_element.get_text() if code_element else ''
    print(text)

// *[ @ id = "supports-class-probabilities"] / ul[2]
// *[ @ id = "supports-class-probabilities"] / ul[2] / li[1]
# hyperparams.append()
print(hyperparams_element)