from django.shortcuts import render

from .models import Item

# Create your views here.
from django.http import HttpResponse
from django.template import loader

from io import StringIO
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def get_gene_prevalence_graph(species_name):
    df = pd.read_csv('./database.csv')

    subdf = df[df['species'] == species_name]

    fig = plt.figure(figsize=(10, 6))
    genes = [*subdf['gene_symbol'].value_counts().head(15).index]
    subdf = subdf[subdf.isin({'gene_symbol': genes})]
    sns.countplot(y=subdf['gene_symbol'].fillna('intergenic'), order=subdf['gene_symbol'].value_counts().index)
    plt.title(f'Gene prevalence in {species_name}')
    plt.xlabel('Number of SNP')
    plt.ylabel('Gene')

    imgdata = StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)

    return imgdata.getvalue()


def get_enrichment_graph(species_name):
    df = pd.read_csv('./database.csv')

    # Очистка gene_symbol
    def clean_gene_symbol(symbols):
        cleaned = []
        for s in symbols:
            if isinstance(s, str):
                # Удаляем аннотации типа '(downstream)' и нечисловые значения
                s_clean = s.split()[0] if ' ' in s else s
                if s_clean.isalpha() or any(c.isalpha() for c in s_clean):
                    cleaned.append(s_clean)
        return list(set(cleaned))

    # Определение видов и их параметров
    species_map = {
        'Gallus gallus': {'gprofiler_id': 'ggallus'},
        'Bos taurus': {'gprofiler_id': 'btaurus'},
        'Homo sapiens': {'gprofiler_id': 'hsapiens'},
        'Mus musculus': {'gprofiler_id': 'mmusculus'},
        'Sus scrofa': {'gprofiler_id': 'sscrofa'},
        'Canis familiaris': {'gprofiler_id': 'cfamiliaris'},
        'Equus caballus': {'gprofiler_id': 'ecaballus'},
        'Oreochromis niloticus': {'gprofiler_id': 'oniloticus'},
    }

    # Получение уникальных генов и gene_id по видам
    gene_data = df.groupby('species').agg({
        'gene_symbol': lambda x: clean_gene_symbol(x.dropna().unique()),
        'gene_id': lambda x: x.dropna().unique().tolist()
    }).to_dict()

    # Функция для анализа обогащения с g:Profiler
    def perform_enrichment_analysis(genes, gene_ids, organism_id, species):
        gp = GProfiler(return_dataframe=True)
        try:
            # Пробуем gene_symbol
            print(f"Попытка анализа обогащения с gene_symbol для {species}: {genes[:5]}...")
            profile = gp.profile(
                organism=organism_id,
                query=genes,
                sources=['GO:BP', 'GO:MF', 'GO:CC'],
                significance_threshold_method='fdr',
                user_threshold=0.05,
                no_evidences=False
            )
            if profile.empty:
                print(f"g:Profiler не нашел аннотации для gene_symbol, пробуем gene_id: {gene_ids[:5]}...")
                profile = gp.profile(
                    organism=organism_id,
                    query=gene_ids,
                    sources=['GO:BP', 'GO:MF', 'GO:CC'],
                    significance_threshold_method='fdr',
                    user_threshold=0.05,
                    no_evidences=False
                )
            if profile.empty:
                print(f"g:Profiler не нашел обогащенных терминов для {species}.")
                return pd.DataFrame()

            # Проверка доступных столбцов
            print(f"Столбцы в ответе g:Profiler: {profile.columns.tolist()}")

            # Фильтрация GO-терминов
            go_terms = profile[profile['source'].str.contains('GO:')][['source', 'name', 'p_value', 'intersection_size']]
            go_terms.rename(columns={'name': 'go_term', 'source': 'go_category'}, inplace=True)
            return go_terms
        except Exception as e:
            print(f"Ошибка g:Profiler для {species}: {e}")
            return pd.DataFrame()

    genes = gene_data['gene_symbol'].get(species_name, [])
    gene_ids = gene_data['gene_id'].get(species_name, [])
    gprofiler_id = species_map[species_name]['gprofiler_id']

    # Анализ обогащения
    enrichment_df = perform_enrichment_analysis(genes, gene_ids, gprofiler_id, species_name)

    if not enrichment_df.empty:
        # Визуализация топ-10 GO-терминов
        fig = plt.figure(figsize=(10, 6))
        sns.barplot(
            y=enrichment_df['go_term'][:10],
            x=-enrichment_df['p_value'][:10].apply(lambda x: -np.log10(x)),
            hue=enrichment_df['go_category'][:10]
        )
        plt.title(f'Top 10 GO terms for {species_name} (-log10 p-value)')
        plt.xlabel('-log10(p-value)')
        plt.ylabel('GO-term')

        plt.legend(title='GO Category')

        imgdata = StringIO()
        fig.savefig(imgdata, format='svg')
        imgdata.seek(0)

        return imgdata.getvalue()

    else:
        print(f"Обогащенные термины для {species_name} не найдены.")


def get_gene_by_chromosome_graph(species_name):
    df = pd.read_csv('./database.csv')

    subdf = df[df['species'] == species_name]

    # Распределение по хромосомам для кур
    fig = plt.figure(figsize=(10, 6))
    sns.countplot(x=subdf['chrom'], hue=subdf['gene_symbol'].fillna('intergenic'))
    plt.title(f'Distribution of genes across chromosomes in {species_name}')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of SNP')

    imgdata = StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)

    return imgdata.getvalue()


def index(request):
    species_list = [*Item.objects.order_by().values_list('species', flat=True).distinct()]
    gene_list = [*Item.objects.order_by('gene_symbol').values_list('gene_symbol', flat=True).distinct()]
    chrom_list = sorted(set([d.get('chrom', '').strip() for d in Item.objects.values_list('data', flat=True)]))
    chrom_list = [v for v in chrom_list if len(v) < 3]

    template = loader.get_template("species/index.html")
    context = {'species_list': species_list, 'gene_list': gene_list, 'chrom_list': chrom_list}

    return HttpResponse(template.render(context, request))


def details(request):
    template = loader.get_template("species/details.html")

    species_name = request.GET.get('species_name', '')
    gene_name = request.GET.get('gene_name', '')

    start_pos = request.GET.get('start_pos', '')
    end_pos = request.GET.get('end_pos', '')

    only_metrics = request.GET.get('only_metrics', '') == 'True'

    qs = Item.objects
    if species_name:
        qs = qs.filter(species=species_name)
    if gene_name:
        qs = qs.filter(gene_symbol=gene_name)

    l = [*qs.all()[:15]]

    if start_pos:
        l = [v for v in l if v.start and start_pos < v.start]

    if end_pos:
        l = [v for v in l if v.start and v.end < end_pos]

    field_names = ['gene_symbol', 'gene_id', 'start', 'end', 'snp_id', 'snp_pos', 'article_title', 'article_url',]
    key_names = [*set(sum([[*ins.data.keys()] for ins in l], []))]
    col_names = field_names + key_names  # output = "<br>".join([q.species + "|" + q.article_title for q in l])
    metric_col_names = ['DDAF', 'EHH', 'Fst', 'H12', 'H2', 'Meta', 'Pi', 'Rsb', 'SHp', 'TajimaD', 'XP', 'XP', 'di', 'iHS', 'metric', 'nSL', 'omega']
    metric_col_names = [
        col_name for col_name in metric_col_names
        if any([ins.data.get(col_name, '') for ins in l])
    ]
    if only_metrics:
        col_names = [col_name for col_name in col_names if col_name in ['gene_symbol', 'gene_id', *metric_col_names]]
    table = [
        [str(getattr(ins, k, '')) for k in field_names if k in col_names]
        + [str(ins.data[k] or '') for k in key_names if k in col_names and not k in field_names]
        for ins in l
    ]
    if table:
        table = [col_names]  + table

    try:
        graph1 = get_gene_prevalence_graph(species_name)
    except Exception as e:
        print('Error: get_gene_prevalence_graph', repr(e))
        graph1 = None

    try:
        graph2 = get_enrichment_graph(species_name)
    except Exception as e:
        print('Error: get_enrichment_graph', repr(e))
        graph2 = None

    try:
        graph3 = get_gene_by_chromosome_graph(species_name)
    except Exception as e:
        print('Error: get_gene_by_chromosome_graph', repr(e))
        graph3 = None

    context = {
        'species_name': species_name,
        "table": table,
        'graph2': graph1,
        'graph3': graph2,
        'graph4': graph3,
        'only_metrics': only_metrics
    }

    return HttpResponse(template.render(context, request))


def items(request):
    template = loader.get_template("species/items.html")

    species_name = request.GET.get('species_name', '')
    gene_name = request.GET.get('gene_name', '')
    chrom_name = request.GET.get('chrom_name', '')

    start_pos = request.GET.get('start_pos', '')
    end_pos = request.GET.get('end_pos', '')

    qs = Item.objects
    if species_name:
        qs = qs.filter(species=species_name)
    if gene_name:
        qs = qs.filter(gene_symbol=gene_name)

    l = [*qs.all()]

    if start_pos:
        l = [v for v in l if v.start and start_pos < v.start]

    if end_pos:
        l = [v for v in l if v.start and v.end < end_pos]

    if chrom_name:
        l = [v for v in l if v.data.get('chrom', '').strip() == chrom_name]

    l = l[:95]

    field_names = ['species', 'gene_symbol', 'gene_id', 'start', 'end', 'snp_id', 'snp_pos', 'article_title', 'article_url',]
    key_names = [*set(sum([[*ins.data.keys()] for ins in l], []))]
    col_names = field_names + key_names
    table = [
        [str(getattr(ins, k, '')) for k in field_names if k in col_names]
        +
        [str(ins.data[k] or '') for k in key_names if k in col_names]
        for ins in l
    ]
    if table:
        table = [col_names]  + table

    graph = None

    context = {"table": table, 'graph': graph}

    return HttpResponse(template.render(context, request))
