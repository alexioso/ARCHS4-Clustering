def purity_to_excel(cluster_metadata_metrics, outfile):
    workbook = Excel.Workbook("tables/" +outfile +".xlsx")
    worksheet = workbook.add_worksheet()
    cluster_format = workbook.add_format({'bold': True, 'align':'center'})
    center = workbook.add_format({'align': 'center'})

    # Add a format. Light red fill with dark red text.
    format1 = workbook.add_format({'bg_color': '#FFC7CE',
                                   'font_color': '#9C0006'})

    # Add a format. Green fill with dark green text.
    format2 = workbook.add_format({'bg_color': '#C6EFCE',
                                   'font_color': '#006100'})

    worksheet.write('A1', "Cluster ID",cluster_format)
    worksheet.write('B1', "Cluster Size",center)
    worksheet.write('C1', "Dominant Molecule",center)
    worksheet.write('D1', "Molecule Purity",center)
    worksheet.write('E1', "Dominant Series",center)
    worksheet.write('F1', "Series Purity",center)
    worksheet.write('G1', "Dominant Label",center)
    worksheet.write('H1', "Label Purity",center)

    worksheet.set_column(1, 7, 15)

    num_keys = len(list(cluster_metadata_metrics.keys()))

    for i in np.sort(list(cluster_metadata_metrics.keys())):
        worksheet.write(i,0,i,center)
        worksheet.write(i,1,cluster_metadata_metrics[i]["Size"],center)
        worksheet.write(i,2,list(cluster_metadata_metrics[i]["Molecule Purity"].keys())[0],center)
        worksheet.write(i,3,list(cluster_metadata_metrics[i]["Molecule Purity"].values())[0],center)
        worksheet.write(i,4,list(cluster_metadata_metrics[i]["Series Purity"].keys())[0],center)
        worksheet.write(i,5,list(cluster_metadata_metrics[i]["Series Purity"].values())[0],center)
        worksheet.write(i,6,list(cluster_metadata_metrics[i]["Label Purity"].keys())[0],center)
        worksheet.write(i,7,list(cluster_metadata_metrics[i]["Label Purity"].values())[0],center)

    # Write a conditional format over a range.
    worksheet.conditional_format('H2:H'+str(num_keys+1), {'type': 'cell',
                                             'criteria': '>=',
                                             'value': 0.8,
                                             'format': format2})

    # Write another conditional format over the same range.
    worksheet.conditional_format('H2:H'+str(num_keys+1), {'type': 'cell',
                                             'criteria': '<',
                                             'value': 0.8,
                                             'format': format1})



    workbook.close()