import os

def generate_html_report(output_dir, results):
    html_path = os.path.join(output_dir, "combination_report.html")
    with open(html_path, 'w') as f:
        f.write("<html><head><title>FASTQ Combiner Report</title>")
        f.write('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">')
        f.write("</head><body class='container mt-4'>")
        f.write("<h1 class='mb-4'>FASTQ Combiner Report</h1>")
        f.write("<table class='table table-striped'><thead><tr>")
        f.write("<th>Sample</th><th>R1 Count</th><th>R2 Count</th><th>Status</th>")
        f.write("</tr></thead><tbody>")
        for res in results:
            f.write(f"<tr><td>{res['sample']}</td><td>{res['r1_count']}</td><td>{res['r2_count']}</td><td>{res['status']}</td></tr>")
        f.write("</tbody></table></body></html>")
    print(f"HTML report written to {html_path}")
