import requests
import sys
from bs4 import BeautifulSoup, NavigableString
from bs4.diagnose import diagnose
import traceback
import collections
import time

collections.Callable = collections.abc.Callable


class Record:
    def __init__(self):
        self.data = []


def searchSpecies(record, soup):
    infoboxbar = soup.find_all(class_="infobox_key")

    # Default, to be overwritten
    record.species = "NULL"

    for infoSection in infoboxbar:
        if "Species: " in infoSection.text:
            record.species = infoSection.text.split(" ")[1] + " " + infoSection.text.split(" ")[2]


def searchTaxonomy(record, soup):
    taxonomyPlace = soup.find(class_="id_1 section expandsection taxonomy-table")

    # Defaults, to be overwritten
    record.taxDomain = "NULL"
    record.taxPhylum = "NULL"
    record.taxClass = "NULL"
    record.taxOrder = "NULL"
    record.taxFamily = "NULL"
    record.taxGenus = "NULL"

    try:
        for section in taxonomyPlace.find_all("tr"):
            sectionLevel = section.find(class_="bold_valigntop width180_valigntop")
            if "Domain" in sectionLevel:
                record.taxDomain = section.find(class_="valigntop paddingright").text
                # print("Domain: " + record.taxDomain)
            if "Phylum" in sectionLevel:
                record.taxPhylum = section.find(class_="valigntop paddingright").text
                # print("Phylum: " + record.taxPhylum)
            if "Class" in sectionLevel:
                record.taxClass = section.find(class_="valigntop paddingright").text
                # print("Class: " + record.taxClass)
            if "Order" in sectionLevel:
                record.taxOrder = section.find(class_="valigntop paddingright").text
                # print("Order: " + record.taxOrder)
            if "Family" in sectionLevel:
                record.taxFamily = section.find(class_="valigntop paddingright").text
                # print("Family: " + record.taxFamily)
            if "Genus" in sectionLevel:
                record.taxGenus = section.find(class_="valigntop paddingright").text
                # print("Genus: " + record.taxGenus)
    except:
        print("No taxonomy table found!")
        pass


def searchPathogenicity(record, soup):
    sections = soup.find_all(class_="resultdetail")

    # Default, to be overwritten
    record.pathogenicity = "False"

    for section in sections:  # Uses tooltip text of section to determine it is that section
        if "Information on possible application of the strain and its possible interaction with e.g. potential hosts" in section.text:
            try:
                for line in section.find_all(class_="mark-hover"):
                    if "Pathogenicity" in line.txt and "yes" in line.find(class_="valigntop"):
                        record.pathogenicity = "True"
                        break  # May find multiple entries; will use the first
            except:
                print("No pathogenicity information found.")
                pass


def searchAerobicity(record, soup):
    sections = soup.find_all(class_="resultdetail")

    # Default, to be overwritten
    record.aerobicity = "NULL"

    for section in sections:  # Uses tooltip text of section to determine it is that section
        if "Information on physiology and metabolism" in section.text:
            try:
                aerobicityLine = section.find(class_="width180_valigntop bold_valigntop", text="Oxygen tolerance")
                if aerobicityLine:
                    record.aerobicity = aerobicityLine.findNext("td").text
            except:
                print("No aerobicity information found.")
                pass


def searchGC_Content(record, soup):
    sections = soup.find_all(class_="resultdetail")

    # Default, to be overwritten
    record.gc_content = "NULL"

    for section in sections:  # Uses tooltip text of section to determine it is that section
        gcSection = section.find(class_="bold_valigntop width180_valigntop", text="GC-content")
        if gcSection:
            try:
                afterElement = gcSection.findNext("td")
                # Converts from percent to decimal
                record.gc_content = str(round(float(afterElement.text.split("mol")[0].strip()) / 100, 4))
                break  # May find multiple entries; will use the first
            except:
                print("No GC Content information found.")
                pass


def searchTagsAndSource(record, soup):
    sections = soup.find_all(class_="resultdetail")

    # Defaults, to be overwritten
    record.tags1 = []
    record.tags2 = []
    record.source = "NULL"

    for section in sections:  # Uses tooltip text of section to determine it is that section
        if "Information on isolation source, the sampling and environmental conditions" in section.text:  # sic
            sourceLine = section.find(class_="bold_valigntop width180_valigntop", text="Sample type/isolated from")
            if sourceLine:
                record.source = sourceLine.findNext("td").text
            try:
                for line in section.find_all(class_="mark-hover"):
                    if "Isolation sources categories" in line.find(class_="bold_valigntop width180_valigntop"):
                        for row in line.find(class_="detail-isol-categories").find_all("tr"):
                            colNum = 1
                            for col in row.select('td[class*="isol-color"]'):
                                if colNum == 1:
                                    record.tags1.append(col.text)
                                elif colNum == 2:
                                    record.tags2.append(col.text)
                                else:
                                    break
                                colNum += 1
            except:
                print("No tags/source found.")
                pass


def populateOneRecord(index):
    url = "https://bacdive.dsmz.de/strain/" + str(index)

    page = requests.get(url)
    if page.status_code == 200:
        print("Successfully read page.")
        record = Record()
        record.isscrape = True
        record.url = url
        record.id = str(index)

        soup = BeautifulSoup(page.text, 'html.parser')
        # print(soup.prettify())

        goodsoup = soup  # Seeming to miss content?
        # goodsoup = soup.find(id="middle_content", recursive=True)
        # goodsoup = soup.find(id="middle_content")

        # print(goodsoup.prettify(), file=sys.stdout)

        searchSpecies(record, goodsoup)
        searchTaxonomy(record, goodsoup)
        # Living state to be determined based on other information.
        searchPathogenicity(record, goodsoup)
        searchAerobicity(record, goodsoup)
        searchGC_Content(record, goodsoup)
        searchTagsAndSource(record, goodsoup)

        if record.pathogenicity == "True" or "Host" in record.tags1:
            record.livingState = "Not Free Living"
        else:
            record.livingState = "Free Living"

            print(record.id + "\t" +
                  record.species + "\t" +
                  record.taxDomain + "\t" +
                  record.taxPhylum + "\t" +
                  record.taxClass + "\t" +
                  record.taxOrder + "\t" +
                  record.taxFamily + "\t" +
                  record.taxGenus + "\t" +
                  record.livingState + "\t" +
                  record.pathogenicity + "\t" +
                  record.aerobicity + "\t" +
                  record.gc_content + "\t" +
                  ", ".join([str(tag) for tag in record.tags1]) + "\t" +
                  ", ".join([str(tag) for tag in record.tags2]) + "\t" +
                  record.source + "\t" +
                  record.url + "\n")

    else:
        print("Status from page " + url + ": " + str(page.status_code))
        raise Exception

    return record


def gatherBacDive():
    print("Gathering BacDive data! (This will take a few minutes.)")


    indices = []
    with open("indexList.txt","r") as indexFile:
        for line in indexFile.readlines():
            indices.append(line.strip())

    i = 0
    with open('bacDiveScrapeResults.txt', 'a+') as outFile:
        for bacDiveIndex in indices:
            try:
                if i > 100:
                    break
                i += 1
                record = populateOneRecord(bacDiveIndex)
                time.sleep(3)

                outFile.write(record.id + "\t" +
                              record.species + "\t" +
                              record.taxDomain + "\t" +
                              record.taxPhylum + "\t" +
                              record.taxClass + "\t" +
                              record.taxOrder + "\t" +
                              record.taxFamily + "\t" +
                              record.taxGenus + "\t" +
                              record.livingState + "\t" +
                              record.pathogenicity + "\t" +
                              record.aerobicity + "\t" +
                              record.gc_content + "\t" +
                              ", ".join([str(tag) for tag in record.tags1]) + "\t" +
                              ", ".join([str(tag) for tag in record.tags2]) + "\t" +
                              record.source + "\t" +
                              record.url + "\n")

            except Exception as e:
                stack = traceback.format_exc()
                print("Sadness: " + str(e) + " (" + str(stack) + ")", file=sys.stderr)

        print("Ended BacDive search.")



gatherBacDive()
