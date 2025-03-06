#!/usr/bin/env python3

import json
import argparse
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

def create_table_data(data, title):
    """Convert dictionary data to table format with title"""
    table_data = [[Paragraph(f"<b>{title}</b>", styles['Heading2']), '']]  # Title row
    table_data.extend([[key, value] for key, value in data.items()])
    return table_data

def create_styled_table(table_data):
    """Create a styled table from data"""
    # Create table with data
    table = Table(table_data, colWidths=[4*inch, 2.5*inch])
    
    # Define style
    style = TableStyle([
        # Title row style
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('SPAN', (0, 0), (-1, 0)),  # Span the title across all columns
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        
        # Content rows style
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
        ('ALIGN', (0, 1), (0, -1), 'LEFT'),  # Left align keys
        ('ALIGN', (1, 1), (1, -1), 'RIGHT'),  # Right align values
        
        # Grid style
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 14),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ])
    table.setStyle(style)
    return table

def create_report(json_data, output_file, image_paths=None):
    """Create PDF report from JSON data and optional images"""
    # Create the document
    doc = SimpleDocTemplate(
        output_file,
        pagesize=letter,
        rightMargin=72,
        leftMargin=72,
        topMargin=72,
        bottomMargin=72
    )
    
    # Initialize story (content list)
    story = []
    
    # Add title
    title = Paragraph("Genome Analysis Report", styles['Title'])
    story.append(title)
    story.append(Spacer(1, 30))
    
    # Parse JSON data
    data = json_data
    
    # Create tables for each section
    for section, content in data.items():
        # Create and add table
        table_data = create_table_data(content, section)
        table = create_styled_table(table_data)
        story.append(table)
        story.append(Spacer(1, 30))
    
    # Add images if provided
    if image_paths:
        for img_path in image_paths:
            try:
                img = Image(img_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 20))
            except Exception as e:
                print(f"Error adding image {img_path}: {e}")
    
    # Build the document
    doc.build(story)

# Initialize styles
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name='CustomTitle',
    parent=styles['Title'],
    fontSize=24,
    spaceAfter=30
))

def main():
    parser = argparse.ArgumentParser(description='Convert JSON to PDF report')
    parser.add_argument('-j', '--json', required=True, help='Path to input JSON file')
    parser.add_argument('-i', '--images', help='Text file with image paths')
    parser.add_argument('-o', '--output', default='AnnoAudit_Report.pdf', help='Output PDF file path')
    
    args = parser.parse_args()

    with open(args.json) as f:
        json_data = json.load(f)
    
    image_paths = []
    if args.images:
        with open(args.images, 'r') as f:
            image_paths = [line.strip() for line in f if line.strip()]
    
    # Create the report
    create_report(json_data, args.output, image_paths)
    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    main()