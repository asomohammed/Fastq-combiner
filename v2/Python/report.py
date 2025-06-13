import os

def generate_html_report(output_dir, results):
    html_path = os.path.join(output_dir, "combination_report.html")
    with open(html_path, 'w') as f:
        f.write("<html><head><title>FASTQ Combiner Report</title>")
        f.write('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">')
        f.write("</head><body class='container mt-4'>")
        f.write("<h1 class='mb-4'>FASTQ Combiner Report</h1>")
        f.write("<table class='table table-striped table-hover'><thead><tr>")
        f.write("<th>Sample</th>")
        f.write("<th>R1 Count</th>")
        f.write("<th>R2 Count</th>")
        f.write("<th>Status</th>")
        f.write("<th>Combined R1 Output</th>")
        f.write("<th>Combined R2 Output</th>")
        f.write("<th>MD5 R1</th>")
        f.write("<th>MD5 R2</th>")
        f.write("</tr></thead><tbody>")

        for res in results:
            f.write("<tr>")
            f.write(f"<td>{res['sample']}</td>")
            f.write(f"<td>{res['r1_count']}</td>")
            f.write(f"<td>{res['r2_count']}</td>")
            f.write(f"<td>{res['status']}</td>")
            f.write(f"<td>{os.path.basename(res['combined_r1_output']) if res['combined_r1_output'] else ''}</td>")
            f.write(f"<td>{os.path.basename(res['combined_r2_output']) if res['combined_r2_output'] else ''}</td>")
            f.write(f"<td>{res['md5_r1']}</td>")
            f.write(f"<td>{res['md5_r2']}</td>")
            f.write("</tr>")

        f.write("</tbody></table></body></html>")

    print(f"HTML report written to {html_path}")
