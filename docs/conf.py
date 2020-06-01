from recommonmark.parser import CommonMarkParser

source_parsers = {
    '.md': CommonMarkParser,
}

extensions = [
    'sphinx_markdown_tables',
]
source_suffix = ['.rst', '.md']

master_doc = 'index'
project = u'iCOMIC_doc'
