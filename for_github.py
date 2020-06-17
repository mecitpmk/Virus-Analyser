from reportlab.lib.units import cm
from Bio import SeqIO
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio.Graphics import BasicChromosome
from tkinter import *
import requests
from bs4 import  BeautifulSoup
import wget,pickle,dbm,os
from PIL import ImageTk
from PIL import Image,ImageDraw
from pdf2image import convert_from_path



class GUI(Frame):
    def __init__(self,parent):
        self.parent=parent
        Frame.__init__(self,parent)
        self.grid()
        self.initGUI()
        self.need_Fetch=True
        if "biology_informations.db.dir" in os.listdir():
            print("Database Founded.")
            self.error_label.config(text="Database Founded",bg="yellow",fg="black")
            self.original_settings(self.error_label)
            self.database_download()
            self.need_Fetch=False
            self.append_to_listbox()
    def initGUI(self):
        self.virus_and_files={}
        a_frame=Frame(self,relief=GROOVE,borderwidth=2)
        a_frame.grid(row=0,column=0,padx=10,pady=10,columnspan=3)
        self.error_label=Label(self,text="",anchor=E,fg="white")
        self.error_label.grid(row=0,column=2,padx=10,pady=10,sticky=E)
        Label(a_frame,text="Click And Collect Viruses Links!").grid(row=0,column=0,padx=10,pady=10,sticky=E)
        self.collect_=Button(a_frame,text="Collect",command=self.capture_datas)
        self.collect_.grid(row=0,column=1,padx=10,pady=10,ipadx=20)
        
        
        self.cvFrame=Frame(self,relief=GROOVE,borderwidth=4)
        self.cvFrame.grid(row=1,column=0,columnspan=3,padx=10,pady=10)
        self.canvas=Canvas(self.cvFrame)
        self.canvas.pack(ipadx=380,ipady=75)


        self.canVeritcalScr=Scrollbar(self.canvas,orient="vertical",command=self.canvas.yview)
        self.canVHorizontalScr=Scrollbar(self.canvas,orient="horizontal",command=self.canvas.xview)
        self.canVHorizontalScr.pack(side="bottom",fill="x")
        self.canVeritcalScr.pack(side="right",fill="y")
        self.canvas.configure(xscrollcommand=self.canVHorizontalScr.set,yscrollcommand=self.canVeritcalScr.set)
        
        Label(self,text="Download Progess Bar").grid(row=2,column=0,columnspan=3,padx=10,pady=10,sticky=S)
        
        self.progress_bar=Canvas(self,bg="grey",width=390,height=20,borderwidth=3,relief="sunken")
        self.progress_bar.grid(row=3,column=0,padx=10,pady=10,columnspan=4)
        
        self.lstboxFrame=Frame(self,relief=GROOVE,borderwidth=2)
        self.lstboxFrame.grid(row=4,column=0,padx=10,pady=10)

        self.all_viruses_listbox=Listbox(self.lstboxFrame,width=30,height=15)
        self.all_viruses_listbox.grid(row=0,column=0,padx=10,pady=10)
        self.all_viruses_listbox.bind("<<ListboxSelect>>",self.clicked_viruses)

        scrool3=Scrollbar(self.lstboxFrame)
        scrool3.grid(row=0,column=0,sticky=E,ipady=97,padx=10,pady=10)
        self.all_viruses_listbox.configure(yscrollcommand=scrool3.set)
        scrool3.config(command=self.all_viruses_listbox.yview)
        

        self.download_viruses_genbnk=Button(self.lstboxFrame,text="Select and\nDownload Virus\nGenbank",command=self.download_virus_bank)
        self.download_viruses_genbnk.grid(row=0,column=1,padx=10,pady=10)

        self.buttons_frames=Frame(self,relief=GROOVE,borderwidth=4)
        self.buttons_frames.grid(row=4,column=1,columnspan=2,padx=10,pady=10)


        Label(self.buttons_frames,text="Functions Buttons").grid(row=0,column=0,columnspan=2)

        self.check_features=Button(self.buttons_frames,text="Check Virus Features",command=self.check_ftrs)
        self.check_features.grid(row=1,column=0,padx=25,pady=25)

        self.check_chromosomes=Button(self.buttons_frames,text="Check Chromosomes",command=self.chromosome_analysis)
        self.check_chromosomes.grid(row=1,column=1,padx=25,pady=25)

        self.check_gene_diagram=Button(self.buttons_frames,text="Check Gene Diagram",command=self.gene_diagram_analysis)
        self.check_gene_diagram.grid(row=2,column=0,padx=25,pady=25)

        self.show_selected_features=Button(self.buttons_frames,text="Show Chromosomes Inside",command=self.check_specified_chromosomes)
        self.show_selected_features.grid(row=2,column=1,padx=25,pady=25)
        self.progress=0
        # self.rects=[]
    def database_download(self):

        db=dbm.open("biology_informations.db","c")
        self.all_virus_dictionary=pickle.loads(db["LINKS"])
        try:
            self.virus_and_files=pickle.loads(db["FILES"])
        except :
            pass
        db.close()

    def upload_database(self,vir_link=False,vir_files=False):

        db=dbm.open("biology_informations.db","c")
        if vir_files:
            db["FILES"]=pickle.dumps(self.virus_and_files)
        elif vir_link:
            db["LINKS"]=pickle.dumps(self.all_virus_dictionary)
        db.close()
    def original_settings(self,label):
        self.after(2000,lambda : label.configure(text="",bg="SystemButtonFace"))
    def capture_datas(self):
        if self.need_Fetch:
            req=requests.get("https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/")
            soup=BeautifulSoup(req.content,'html.parser')
            all_a=soup.find_all('a')
            self.viruses={}
            vir_url="https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/"
            for a in all_a:
                if a.text != "DBV/" and a.text != "Parent Directory" and a.text != "FamilyPhylogeneticTree/" and a.text != "AlignFastaV/" and a.text[0:4] != "all." \
                    and a.text != "Viruses_RefSeq_and_neighbors_genome_data.tab":
                    self.viruses[a.text]=f'{vir_url}{a.text}'
                    self.update()
            self.clear_datas()
        else:
            print("Already Fetched.")
    def bar_custom(self,current,total,width=80):
        print("Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total))
        try:
            self.progress=(current/total)*400
        except:
            pass
        else:
            self.test=self.progress_bar.create_rectangle(0,0,self.progress,30,fill="green")
            self.update()
    def download_virus_bank(self):
        print(self.selected_link)
        req=requests.get(self.selected_link)
        soup=BeautifulSoup(req.content,'html.parser')
        will_download=soup.find_all('a')
        downloaded=False
        if self.selected_viruses not in self.virus_and_files:
            self.virus_and_files[self.selected_viruses]=[]
        for links in will_download:
            splitted=links.text.split(".")
            for text in splitted:
                if text == 'gbk':
                    if links.text not in self.virus_and_files[self.selected_viruses]:
                        self.progress_bar.delete("all")
                        self.update()
                        link=wget.download(f'{self.selected_link}{links.text}',bar=self.bar_custom)
                        # print(link)
                        self.progress=0
                        self.virus_and_files[self.selected_viruses].append(link)
                        print(self.virus_and_files[self.selected_viruses])
                        self.error_label.config(text="Download Complete!",bg="green",fg="white")
                        self.original_settings(self.error_label)
                        downloaded=True
                    else:
                        self.progress_bar.delete("all")
                        self.progress_bar.create_rectangle(0,0,400,30,fill="yellow")
            if downloaded:
                self.upload_database(vir_files=True)
    def clear_datas(self):
        self.all_virus_dictionary={}
        for virus,virus_link in self.viruses.items():
            empty_str=""
            virus_str=virus
            virus_str=virus_str.strip("/").split("_")
            for names in virus_str:
                if names[0:3] != "uid":
                    empty_str=f'{empty_str}{names} '
            self.all_virus_dictionary.setdefault(empty_str.capitalize(),self.viruses[virus])
        self.append_to_listbox()
        self.upload_database(vir_link=True)
        
    def append_to_listbox(self):
        for virus_names in self.all_virus_dictionary.keys():
            self.all_viruses_listbox.insert(END,virus_names)
    
    def clicked_viruses(self,event):
        capture_selected=event.widget
        selected=capture_selected.curselection()
        try:
            self.selected_viruses=self.all_viruses_listbox.get(selected)
            self.selected_link=self.all_virus_dictionary[self.selected_viruses]
        except :
            pass
    def check_ftrs(self):
        if hasattr(self,"selected_viruses"):
            if self.selected_viruses in self.virus_and_files:
                print("Okey .")
                self.progress_bar.create_rectangle(0,0,400,30,fill="yellow")
                self.show_ftrs(virus_name=self.selected_viruses)
                self.top_lv=Toplevel(self)
                self.top_lv.title(string=self.selected_viruses)
                Label(self.top_lv,text=f'Virus Name: {self.selected_viruses}',bg="yellow",font=("Helvetica",12,"bold")).grid(row=0,column=0,columnspan=4,padx=10,pady=10)
                Label(self.top_lv,text="All Genes Include:").grid(row=1,column=0,padx=10,pady=10)
                genes_lb=Listbox(self.top_lv)
                genes_lb.grid(row=1,column=1,padx=10,pady=10)

                Label(self.top_lv,text="Locus Tags: ").grid(row=2,column=0,padx=10,pady=10)
                locus_lb=Listbox(self.top_lv)
                locus_lb.grid(row=2,column=1,padx=10,pady=10)


                scrool=Scrollbar(self.top_lv)
                scrool.grid(row=1,column=1,sticky=E,ipady=55,padx=10,pady=10)
                genes_lb.configure(yscrollcommand=scrool.set)
                scrool.config(command=genes_lb.yview)


                scrool2=Scrollbar(self.top_lv)
                scrool2.grid(row=2,column=1,sticky=E,ipady=55,padx=10,pady=10)
                locus_lb.configure(yscrollcommand=scrool2.set)
                scrool2.config(command=locus_lb.yview)

                for g in self.genes_list: genes_lb.insert(END,g)
                for l in self.locus_tags_list: locus_lb.insert(END,l)

                frame_text=Frame(self.top_lv,relief=GROOVE,borderwidth=3)
                frame_text.grid(row=1,column=2,rowspan=2,padx=10,pady=10)

                frame_descript=Frame(self.top_lv,relief=GROOVE,borderwidth=3)
                frame_descript.grid(row=1,column=3,rowspan=2,padx=10,pady=10)

                c_label=Label(frame_text)
                c_label.grid(row=0,column=0,padx=10,pady=10)
                if len(self.collect_date) > 0:
                    c_label.configure(text=f'Collected Date : {self.collect_date[0]}')
                else:
                    c_label.configure(text=f'Collected Date : Undefined')
                
                c_ctry=Label(frame_text)
                c_ctry.grid(row=1,column=0,padx=10,pady=10)
                if len(self.collect_country) > 0 :
                    c_ctry.configure(text=f"Collected Country : {self.collect_country[0]}")
                else:
                    c_ctry.configure(text=f'Collected Country : Undefined')
                text_mol="Molecul Type : Undefined"
                if len(self.mol_type) >= 1:
                    text_mol="Molecul Type : "
                    for a in range(len(self.mol_type)):
                        text_mol+=f'{a+1} - {self.mol_type[a]}\n'
                Label(frame_text,text=text_mol).grid(row=2,column=0,padx=10,pady=10)
                text_sg="Segment Type : Undefined"
                if len(self.segment) >= 1 :
                    text_sg="Segment Type : \n"
                    for a in range(len(self.segment)):
                        text_sg+=f'{a+1} - {self.segment[a]}\n'
                Label(frame_text,text=text_sg).grid(row=3,column=0,padx=10,pady=10)
                
                self.Seq_ShowButton=Button(frame_text,text="Show Gene Sequence",command=self.show_g_seq)
                self.Seq_ShowButton.grid(row=4,column=0,padx=10,pady=10)
                
                Label(frame_descript,text=f'Description : \n{self.vir_features[self.selected_viruses]["Description"]}').grid(row=0,column=0,padx=10,pady=10)
                str_for_Taxonomy=""
                for tax in self.vir_features[self.selected_viruses]["Taxonomy"]:
                    str_for_Taxonomy+=f'{tax}\n'
                Label(frame_descript,text=f"Taxonomy : \n{str_for_Taxonomy}").grid(row=1,column=0,padx=10,pady=10)

                Label(frame_descript,text=f'Topology: \n{self.cur_topology}').grid(row=2,column=0,padx=10,pady=10)
                Label(frame_descript,text=f'Exact Molecule Type: \n{self.cur_mol_t}').grid(row=3,column=0,padx=10,pady=10)


            else:
                self.progress_bar.delete("all")
                self.progress_bar.create_rectangle(0,0,400,30,fill="red")

                print("Download Files First!")
                self.error_label.config(text="You Should Download Files First! ",bg="red",fg="white")
                self.original_settings(self.error_label)
        else:
            self.error_label.config(text="You Should Select Viruses First!!",bg="red")
            self.original_settings(self.error_label)


    def on_configure(self,event):
        canvas_=event.widget
        # update scrollregion after starting 'mainloop'
        # when all widgets are in canvas
        canvas_.configure(scrollregion=canvas_.bbox('all'))

    def on_c(self,event):
        canvas_width=event.width
        self.canvas_.itemconfig(self.a,width=canvas_width)
    
    def show_g_seq(self):
        genseq_toplevel=Toplevel(self.top_lv)
        genseq_toplevel.title(string=self.selected_viruses)

        self.canvas_=Canvas(genseq_toplevel,height=900,width=570)
        self.canvas_.pack(side=LEFT,fill=BOTH,expand=True)

        sc_bar=Scrollbar(genseq_toplevel,command=self.canvas_.yview,orient="vertical")
        sc_bar.pack(side=RIGHT,fill="y")

        self.canvas_.configure(yscrollcommand=sc_bar.set)

        self.canvas_.bind('<Configure>', self.on_configure)
        
        fram=Frame(self.canvas_)
        
        self.a=self.canvas_.create_window((4,4),window=fram,anchor=NW)



        for files in self.sequence_dict:
            s_dict=self.sequence_dict[files]
            seq=s_dict["Seq"]
            ct=1
            other_ct=0
            texting=""
            number=1
            row_num=0
            all_ct=0
            for a in seq:
                if ct == 1:
                    Label(fram,text=str(number),font=("Helvatica",8,"bold")).grid(row=row_num,column=0)
                    texting+=f'{a}'
                    ct+=1
                    other_ct+=1
                    all_ct+=1
                    self.canvas_.update_idletasks()
                elif ct == 60:
                    texting+=f'{a}'
                    ct=1
                    other_ct=0
                    number+=60
                    all_ct+=1
                    b=Entry(fram,width=520)
                    b.grid(row=row_num,column=1,sticky=W)
                    b.insert(END,texting)
                    b.config(state="readonly")

                    row_num+=1
                    texting=''
                    self.canvas_.update_idletasks()
                    b.update_idletasks()
                elif other_ct == 10:
                    texting+=f'    {a}'
                    other_ct = 1
                    ct+=1
                    all_ct+=1
                    self.canvas_.update_idletasks()
                else:
                    texting+=f'{a}'
                    ct+=1
                    other_ct+=1
                    all_ct+=1
                    self.update_idletasks()
                    if all_ct == len(seq):
                        b=Entry(fram,width=520)
                        b.grid(row=row_num,column=1,sticky=W)
                        b.insert(END,texting)
                        b.config(state="readonly")
                        b.update_idletasks()
            break
    

    def show_ftrs(self,virus_name):
        entries= [(virus_name,self.virus_and_files[virus_name][0])]
        vf=self.virus_and_files[self.selected_viruses]
        entries=[]
        for vfil in vf:
            entries.append((virus_name,vfil))
        self.genes_list=[]
        self.locus_tags_list=[]
        self.collect_date=[]
        self.collect_country=[]
        self.mol_type=[]
        self.segment=[]
        self.sequence_dict={}
        self.vir_features={}
        for vir_name , vir_file  in entries:
            record=SeqIO.read(vir_file,"genbank")
            
            if self.selected_viruses not in self.vir_features:
                self.vir_features.setdefault(self.selected_viruses,{})
                current_virus_dict=self.vir_features[self.selected_viruses]
                current_virus_dict["Taxonomy"]=[]
                current_virus_dict["Description"]=""
                for tax in record.annotations["taxonomy"]:
                    if tax not in current_virus_dict["Taxonomy"]:
                        current_virus_dict["Taxonomy"].append(tax)
                if record.description != current_virus_dict["Description"]:
                    current_virus_dict["Description"]=record.description
                
 
            self.cur_topology=record.annotations.get("topology") 
            self.cur_mol_t=record.annotations.get("molecule_type") 


            if len(entries) >= 1:
                for indx in range(len(entries)):
                    self.sequence_dict[f'{indx+1}st File Sequence']={"Seq":record.seq,"Seq Complement":record.seq.complement()}

            for x in record.features:
                if x.type == "source":
                    if x.qualifiers.get("collection_date") != None:
                        if x.qualifiers.get("collection_date")[0] not in self.collect_date:
                            self.collect_date.append(x.qualifiers.get("collection_date")[0])
                    if x.qualifiers.get("country") != None:
                        if x.qualifiers.get("country")[0] not in self.collect_country:
                            self.collect_country.append(x.qualifiers.get("country")[0]) 
                    if x.qualifiers.get("mol_type") != None:
                        if x.qualifiers.get("mol_type")[0] not in self.mol_type:
                            self.mol_type.append(x.qualifiers.get("mol_type")[0])
                    if x.qualifiers.get("segment") != None:
                        if x.qualifiers.get("segment")[0] not in self.segment:
                            self.segment.append(x.qualifiers.get("segment")[0])

                    break
 
            virus_genes=[feature for feature in record.features if feature.type == "gene"]
            for gene in virus_genes:
                # print(gene)
                for quaf in gene.qualifiers:
                    gen=gene.qualifiers.get("gene")
                    locus=gene.qualifiers.get("locus_tag")
                    if gen  == None:
                        gen=['']
                    elif locus == None:
                        locus=['']
                    for g,l in zip(gen,locus):
                        if g not in self.genes_list and g != '':
                            self.genes_list.append(g)
                        if l not in self.locus_tags_list and l != '':
                            self.locus_tags_list.append(l)
        
    def chromosome_analysis(self):
        if hasattr(self,"selected_viruses"):
            if self.selected_viruses in self.virus_and_files:
                entries=[]
                v_file=self.virus_and_files[self.selected_viruses]
                for f in v_file:
                    entries.append((self.selected_viruses,f))

                # entries = [
                #     (self.selected_viruses, self.virus_and_files[self.selected_viruses]),]

                entry=[]
                for index,(name,filename) in enumerate(entries):
                    record= SeqIO.read(filename, "genbank")
                    entry.append((name,len(record)))

                max_len=max(entry,key=lambda x : x[1])[1]

                telomere_length = 1000 #For illustration

                chr_diagram = BasicChromosome.Organism(output_format="png")
                chr_diagram.page_size = (29.7*cm, 21*cm) #A4 landscape

                for index, (name,leng) in enumerate(entry):
                    # record = SeqIO.read(filename, "fasta")
                    length=leng
                    cur_chromosome = BasicChromosome.Chromosome(name)
                    #Set the scale to the MAXIMUM length plus the two telomeres in bp,
                    #want the same scale used on all five chromosomes so they can be
                    #compared to each other
                    cur_chromosome.scale_num = max_len + 2 * telomere_length

                    #Add an opening telomere
                    start = BasicChromosome.TelomereSegment()
                    start.scale = telomere_length
                    cur_chromosome.add(start)

                    #Add a body - using bp as the scale length here.
                    body = BasicChromosome.ChromosomeSegment()
                    body.scale = length
                    cur_chromosome.add(body)

                    #Add a closing telomere
                    end = BasicChromosome.TelomereSegment(inverted=True)
                    end.scale = telomere_length
                    cur_chromosome.add(end)

                    #This chromosome is done
                    chr_diagram.add(cur_chromosome)

                chr_diagram.draw(f'{self.selected_viruses}.png', self.selected_viruses)
                self.created_png=f'{self.selected_viruses}.png'
                self.show_in_canvas(im_names=self.created_png)
            else:
                self.error_label.config(text="You Should Download Files First!",bg="red")
                self.original_settings(self.error_label)
        else:
            self.error_label.configure(text="You Should Select Viruses First!!",bg="red")
            self.original_settings(self.error_label)
    def show_in_canvas(self,im_names):
        im = Image.open(im_names)
        self.canvas.image = ImageTk.PhotoImage(im)
        self.canvas.create_image(0,0, anchor=NW ,image=self.canvas.image)
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))


    def gene_diagram_analysis(self):
        if hasattr(self,"selected_viruses"):
            if self.selected_viruses in self.virus_and_files:
                gd_diagram = GenomeDiagram.Diagram(self.selected_viruses)
                gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
                gd_feature_set = gd_track_for_features.new_set()

                virfl=self.virus_and_files[self.selected_viruses]
                for vfl in virfl:
                    record=SeqIO.read(vfl,"genbank")

                    for feature in record.features:
                        if feature.type == "gene":
                            #Exclude this feature
                            # continue
                            if len(gd_feature_set) % 2 == 0:
                                color = colors.blue
                            else:
                                color = colors.lightblue
                            gd_feature_set.add_feature(feature,colour=color, label=True,altcolour=color.red)

                gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                                start=0, end=len(record), circle_core=0.72)
                gd_diagram.write(f'{self.selected_viruses}_diagrams.png', "PNG")

                self.diagram_file_names=f'{self.selected_viruses}_diagrams.png'
                self.show_in_canvas(im_names=self.diagram_file_names)
            else:
                self.error_label.config(text="You Should Download Files First!",bg="red")
                self.original_settings(self.error_label)
        else:
            self.error_label.config(text="You Should Select Viruses First!!",bg="red")
            self.original_settings(self.error_label)
    def check_specified_chromosomes(self):
        if hasattr(self,"selected_viruses"):
            if self.selected_viruses in self.virus_and_files:
                v_fil=self.virus_and_files[self.selected_viruses]
                entries=[]
                
                ctr=0
                for v in v_fil:
                    entries.append((f'{self.selected_viruses}{ctr}',v))
                    ctr+=1
                
                max_len=max(entries,key=lambda x : x[1])[1]

                # max_len = 30432563  # Could compute this from the entries dict
                telomere_length = 1000000  # For illustration

                chr_diagram = BasicChromosome.Organism()
                chr_diagram.page_size = (29.7*cm, 21*cm)  # A4 landscape

                for index, (name, filename) in enumerate(entries):
                    record = SeqIO.read(filename, "genbank")
                    max_len=len(record)
                    # for x in record.features:
                    #     print(x.type)
                    length = len(record)
                    features = [f for f in record.features if f.type == "gene"] #tRNA
                    # Record an Artemis style integer color in the feature's qualifiers,
                    # 1 = Black, 2 = Red, 3 = Green, 4 = blue, 5 =cyan, 6 = purple
                    for f in features:
                        f.qualifiers["color"] = [index + 2]

                    cur_chromosome = BasicChromosome.Chromosome(name)
                    # Set the scale to the MAXIMUM length plus the two telomeres in bp,
                    # want the same scale used on all five chromosomes so they can be
                    # compared to each other
                    cur_chromosome.scale_num = max_len + 2 * telomere_length

                    # Add an opening telomere
                    start = BasicChromosome.TelomereSegment()
                    start.scale = telomere_length
                    cur_chromosome.add(start)

                    # Add a body - again using bp as the scale length here.
                    body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
                    body.scale = length
                    cur_chromosome.add(body)

                    # Add a closing telomere
                    end = BasicChromosome.TelomereSegment(inverted=True)
                    end.scale = telomere_length
                    cur_chromosome.add(end)

                    # This chromosome is done
                    chr_diagram.add(cur_chromosome)
                    chr_diagram.draw(f"{self.selected_viruses}chrom.pdf", self.selected_viruses)
                    pages=convert_from_path(f'{self.selected_viruses}chrom.pdf',dpi=100)
                    for page in pages:
                        page.save(f'{self.selected_viruses}chrom.jpg',format="JPEG")
                    self.show_in_canvas(f'{self.selected_viruses}chrom.jpg')
            else:
                self.error_label.config(text="You Should Download Files First!",bg="red")
                self.original_settings(self.error_label)
        else:
            self.error_label.config(text="You Should Select Viruses First!!",bg="red",fg="black")
            self.original_settings(self.error_label)
def main():        
    root=Tk()
    rungui=GUI(root)
    root.title("Analyser Program")
    root.mainloop()
if __name__ == "__main__":
    main()