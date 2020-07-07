#sudo pip install splinter
#sudo pip install pandas


from splinter import Browser
import pandas
import sys
import os
import time

org=sys.argv[1].strip()

browser=Browser(driver_name="chrome",headless=True)
browser.visit("https://www-is.biotoul.fr/search.php")
browser.fill('host', org)
browser.find_by_name('Onsubmit').click()
f=pandas.read_html(browser.find_by_tag('section').html,header=0,index_col=0)[0]

# store IS information??

LLL=f['Name'].tolist()

#print LLL

browser.quit()


date=str(time.strftime("%d%m%y"))
outname="IS_%s_%s.fasta" %(org,date)
outF=open(outname,"w")

for i in LLL:
        i=i.strip()
        address="https://www-is.biotoul.fr/scripts/ficheIS.php?name="+str(i)
        browser=Browser(driver_name="chrome",headless=True)
        browser.visit(address)
        g=str(browser.find_by_id("seq_ident").html).split('entete_propriete">')
        FAM=g[1].split("span>")[1].split("<")[0]
        GRO=g[2].split("span>")[1].split("<")[0]
        seq=str(browser.find_by_tag("div")[11].html).split('seq">')[1].split(" </div")[0].replace("<br>\n","")

        outF.write(">%s\tFamily:%s\tGroup:%s\n%s\n" %(i,FAM,GRO,seq))

        browser.quit()

outF.close()
