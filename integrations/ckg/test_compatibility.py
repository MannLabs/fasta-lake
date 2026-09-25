"""Regression tests for the isolated CKG plotting/network compatibility patch."""
import copy
import sys
import unittest
import os
if os.environ.get("CKG_SOURCE"):
    sys.path.insert(0, os.environ["CKG_SOURCE"])
import networkx as nx
import pandas as pd
from ckg.analytics_core.analytics import analytics
from ckg.analytics_core.viz import viz
from ckg.analytics_core import utils

class NetworkCompatibility(unittest.TestCase):
    def test_edge_callback_removes_a_real_bridge(self):
        graph = nx.path_graph(4)
        self.assertEqual(analytics.most_central_edge(graph), (1, 2))
        components = next(nx.community.girvan_newman(graph, most_valuable_edge=analytics.most_central_edge))
        self.assertEqual({frozenset(c) for c in components}, {frozenset({0,1}),frozenset({2,3})})

    def test_edge_callback_rejects_edgeless_input(self):
        with self.assertRaisesRegex(ValueError, "edgeless"):
            analytics.most_central_edge(nx.empty_graph(2))

    def test_greedy_communities_preserve_disconnected_nodes(self):
        graph = nx.Graph([(1,2),(3,4)])
        result = analytics.get_network_communities(graph, {"communities_algorithm":"greedy_modularity", "values":"width"})
        self.assertEqual(set(result), set(graph))
        self.assertEqual(result[1],result[2])
        self.assertEqual(result[3],result[4])
        self.assertNotEqual(result[1],result[3])

    def test_networkx_table_export_preserves_attributes(self):
        graph = nx.Graph()
        graph.add_node("peptide:P", kind="peptide")
        graph.add_node("protein:P", kind="protein")
        graph.add_edge("peptide:P", "protein:P", width=2, relation="reported membership")
        nodes,edges = viz.network_to_tables(graph,"source","target")
        self.assertEqual(set(nodes["kind"]), {"peptide","protein"})
        self.assertEqual(edges.iloc[0]["relation"], "reported membership")
        self.assertEqual(edges.iloc[0]["width"],2)

    def test_dataframe_network_roundtrip_with_current_pandas(self):
        data=pd.DataFrame({"source":["peptide:A","peptide:B","peptide:B"], "target":["protein:P1","protein:P1","protein:P2"], "width":[1.,1.,1.]})
        result=viz.get_network(data.copy(),"test",{"source":"source","target":"target","values":"width","communities_algorithm":"greedy_modularity","color_weight":False,"limit":None})
        nodes,edges=result["net_tables"]
        self.assertEqual(len(nodes),4)
        self.assertEqual(len(edges),3)
        self.assertIsNotNone(result["app"])

    def test_pyvis_defaults_keep_graph_and_arguments_unchanged(self):
        graph=nx.Graph()
        graph.add_edge("a","b",weight=3, width=2)
        graph.add_node("a",label="Peptide")
        before=copy.deepcopy(graph)
        args={}
        view=viz.get_notebook_network_pyvis(graph,args)
        self.assertEqual(args,{})
        self.assertEqual(dict(graph.nodes(data=True)),dict(before.nodes(data=True)))
        self.assertEqual(dict(((a,b),d) for a,b,d in graph.edges(data=True)),dict(((a,b),d) for a,b,d in before.edges(data=True)))
        self.assertEqual(view.width,"800px")
        self.assertEqual(view.height,"850px")
        self.assertGreater(len(view.html),1000)

    def test_pyvis_css_dimensions_are_not_swapped(self):
        view=viz.get_notebook_network_pyvis(nx.path_graph(2),{"width":"90%","height":"420px"})
        self.assertEqual(view.width,"90%")
        self.assertEqual(view.height,"420px")

    def test_pyvis_current_renderer_retains_html(self):
        view=viz.get_notebook_network_pyvis(nx.path_graph(2),{"width":600,"height":300})
        result=utils.generate_html(view)
        self.assertIsInstance(result,str)
        self.assertEqual(result,view.html)
        self.assertIn("vis.Network",result)

    def test_bar_chart_preserves_counts(self):
        data=pd.DataFrame({"kind":["peptides","proteins"],"count":[3,2]})
        view=viz.get_barplot(data,"counts",{"x":"kind","y":"count"})
        self.assertEqual(list(view.figure["data"][0].y),[3,2])

if __name__ == "__main__":
    unittest.main(verbosity=2)
