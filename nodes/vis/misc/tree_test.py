import json
import re

from more_itertools import chunked

from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle

# def layout(node):
#     if node.is_leaf():
#         N = AttrFace("name", fsize=30)
#         faces.add_face_to_node(N, node, 0, position="aligned")

def get_example_tree():

    # Set dashed blue lines in all leaves
    # nst1 = NodeStyle()
    # nst1["bgcolor"] = "LightSteelBlue"
    # nst2 = NodeStyle()
    # nst2["bgcolor"] = "Moccasin"
    # nst3 = NodeStyle()
    # nst3["bgcolor"] = "DarkSeaGreen"
    # nst4 = NodeStyle()
    # nst4["bgcolor"] = "Khaki"


    t = Tree("((((a1,a2),a3), ((b1,b2),(b3,b4))), ((c1,c2),c3));")
    # t = Tree("((6669.DappuP312785:1.24473,(((7739.JGI126010:4.02e-06,7739.JGI126021:0.168081)0.99985:0.848895,((45351.NEMVEDRAFT_v1g217973-PA:1.49614,(6085.XP_002167371:1.08059,10228.TriadP54105:1.08509)0.457684:0.128729)0.99985:0.321915,(7668.SPU_016633tr:0.00222941,7668.SPU_013365tr:0.093091)0.99985:0.642533)0.888302:0.101463)0.998565:0.12537,((51511.ENSCSAVP00000011400:0.259936,7719.ENSCINP00000035803:0.209259)0.99985:1.762,(7757.ENSPMAP00000006833:0.945887,((7955.ENSDARP00000103018:0.187989,(8049.ENSGMOP00000003903:0.283082,((8128.ENSONIP00000009923:0.172432,(69293.ENSGACP00000005927:0.11289,(99883.ENSTNIP00000008530:0.0637802,31033.ENSTRUP00000029741:0.0851327)0.99985:0.111539)0.750558:0.0140823)0.99985:0.0210439,(8083.ENSXMAP00000007211:0.156672,8090.ENSORLP00000021771:0.242897)0.99985:0.0332475)0.99985:0.0502247)0.99985:0.097906)0.99985:0.204816,(7897.ENSLACP00000021045:0.20103,(8364.ENSXETP00000053091:0.389689,((9258.ENSOANP00000015194:0.652603,((13616.ENSMODP00000032296:0.0669026,(9315.ENSMEUP00000005072:0.0368095,9305.ENSSHAP00000018667:0.0280404)0.997706:0.0140274)0.99985:0.109803,((9361.ENSDNOP00000011178:0.093465,(9371.ENSETEP00000005953:0.276816,(9813.ENSPCAP00000009886:0.0800451,9785.ENSLAFP00000023898:0.0550027)0.99697:0.0294647)0.99985:0.0375967)0.608183:0.00574174,((((30608.ENSMICP00000006810:0.0353003,30611.ENSOGAP00000013460:0.0515798)0.99985:0.034051,(9478.ENSTSYP00000011163:0.0656837,(9483.ENSCJAP00000037059:0.0556713,(9544.ENSMMUP00000001688:0.00821597,(61853.ENSNLEP00000001585:0.00543935,(9601.ENSPPYP00000009641:0.00544861,(9593.ENSGGOP00000014253:0.00241737,(9606.ENSP00000299886:0.00108745,9598.ENSPTRP00000016297:0.00216932)0.993995:0.00108247)0.99985:0.00217395)0.99985:0.00216928)0.99985:0.00380297)0.99985:0.00880979)0.99985:0.0259329)0.993631:0.00725957)0.99985:0.00775695,((((10141.ENSCPOP00000018381:0.0312261,10141.ENSCPOP00000003239:0.0560276)0.99985:0.131426,(10090.ENSMUSP00000018805:0.0326181,10116.ENSRNOP00000003804:0.0330137)0.99985:0.081468)0.995671:0.0105266,(43179.ENSSTOP00000001636:0.0584822,10020.ENSDORP00000002346:0.148367)0.976029:0.00712029)0.99985:0.020967,(9986.ENSOCUP00000010215:0.342444,37347.ENSTBEP00000010812:0.0710806)0.981559:0.0148504)0.99985:0.0111655)0.99985:0.0141806,((9823.ENSSSCP00000018275:0.0497089,(9739.ENSTTRP00000004916:0.0624029,9913.ENSBTAP00000007999:0.0933538)0.99985:0.0153799)0.99985:0.0353429,((9796.ENSECAP00000006507:0.0593149,(9685.ENSFCAP00000004377:0.0895056,(9615.ENSCAFP00000006718:0.0344248,(9646.ENSAMEP00000002816:0.0246918,9669.ENSMPUP00000014182:0.157245)0.979337:0.0106908)0.990824:0.00781435)0.99985:0.0209773)0.99291:0.00461526,(59463.ENSMLUP00000015199:0.175449,132908.ENSPVAP00000006353:0.0636063)0.828177:0.00861023)0.805512:0.00313751)0.99985:0.0138134)0.99985:0.0180555)0.99985:0.0925821)0.99985:0.0504324)0.99985:0.0408666,(28377.ENSACAP00000014302:0.2972,(13735.ENSPSIP00000020553:0.125656,(59729.ENSTGUP00000002997:0.21101,(9103.ENSMGAP00000006649:0.034221,9031.ENSGALP00000007041:0.0301602)0.99985:0.14558)0.99985:0.0950631)0.912514:0.0193477)0.99985:0.0543705)0.99985:0.0558399)0.99985:0.0851946)0.996306:0.0767894)0.99985:0.242292)0.99843:0.214734)0.760404:0.170331)0.99985:0.889243)0.99985:0.309174,((7070.TC009561-PA:1.414,(7029.ACYPI001869-PA:2.67503,((34740.HMEL012959-PA:0.216556,(13037.EHJ66433:0.30242,7091.BGIBMGA012450-TA:0.234676)0.691542:0.0356537)0.99985:1.06833,(7425.NV10250-PA:0.448078,(7460.GB18353-PA:0.265179,12957.ACEP_00009457-PA:0.219522)0.99985:0.223062)0.99985:0.715703)0.950323:0.0960579)0.630765:0.0828761)0.868488:0.104171,(121225.PHUM413450-PA:1.38201,(((43151.ADAR004533-PA:0.286189,7165.AGAP001322-PA:0.216522)0.99985:0.404262,(7176.CPIJ007948-PA:0.310705,7159.AAEL007141-PA:0.261808)0.99985:0.238196)0.99985:0.352684,(7260.FBpp0240708:0.24283,((7222.FBpp0149372:0.191481,7244.FBpp0224481:0.114954)0.99985:0.156407,(7237.FBpp0281462:0.118908,(7217.FBpp0120932:0.132407,(7245.FBpp0269552:0.0320153,7227.FBpp0082093:0.0340969)0.99985:0.0854297)0.99985:0.0800303)0.99985:0.0655239)0.825983:0.0710402)0.99985:1.14395)0.99985:0.708932)0.628915:0.0819834)0.99985:0.309174);")
    # for n in t.traverse():
    #     n.dist = 0

    # n1 = t.get_common_ancestor("a1", "a2", "a3")
    # n1.set_style(nst1)
    # n2 = t.get_common_ancestor("b1", "b2", "b3", "b4")
    # n2.set_style(nst2)
    # n3 = t.get_common_ancestor("c1", "c2", "c3")
    # n3.set_style(nst3)
    # n4 = t.get_common_ancestor("b3", "b4")
    # n4.set_style(nst4)
    ts = TreeStyle()
    # ts.layout_fn = layout
    # ts.show_leaf_name = False

    ts.mode = "c"
    # ts.root_opening_factor = 1
    return t, ts

if __name__ == "__main__":
    t, ts = get_example_tree()
    #t.render("node_background.png", w=400, tree_style=ts)
    ts.force_topology = True
    # t.show(tree_style=ts)
    t.render("nodes/vis/misc/mytree.svg", w=183, units="mm", tree_style=ts)

    x1, y1 = [], []
    x2, y2 = [], []
    x3, y3 = [], []
    with open("nodes/vis/misc/mytree.svg") as f:
        text = f.readlines()[0].replace("</g>", "</g>\n")
        for l in text.split("</g>")[2:]:
            tmp = l.lstrip().rstrip()
            # print(tmp)

            if "polyline" in tmp:
                oldX1, oldY1, oldX2, oldY2 = re.findall('points="(.*?),(.*?)\s(.*?),(.*?)\s"', tmp)[0]
                a, b, c, d, e, f = re.findall('transform="matrix\((.*?),(.*?),(.*?),(.*?),(.*?),(.*?)\)', tmp)[0]
                newX1 = float(a) * float(oldX1) + float(c) * float(oldY1) + float(e)
                newY1 = float(b) * float(oldX1) + float(d) * float(oldY1) + float(f)
                newX2 = float(a) * float(oldX2) + float(c) * float(oldY2) + float(e)
                newY2 = float(b) * float(oldX2) + float(d) * float(oldY2) + float(f)
                # if newX2 - newX1 < 1:
                #     newX2 = newX2 - abs(newX2 - newX1)*0.3
                # else:
                #     newX2 = newX2 + abs(newX2 - newX1)*0.3
                # if newY2 - newY1 < 1:
                #     newY2 = newY2 - abs(newY2 - newY1)*0.3
                # else:
                #     newY2 = newY2 + abs(newY2 - newY1)*0.3
                x1 += [newX1, newX2, None]
                y1 += [newY1, newY2, None]

            if "circle" in tmp:
                oldX, oldY = re.findall('circle cx="(.*?)" cy="(.*?)"', tmp)[0]
                a, b, c, d, e, f = re.findall('transform="matrix\((.*?),(.*?),(.*?),(.*?),(.*?),(.*?)\)', tmp)[0]
                newX = float(a) * float(oldX) + float(c) * float(oldY) + float(e)
                newY = float(b) * float(oldX) + float(d) * float(oldY) + float(f)
                x2 += [newX]
                y2 += [newY]

            if "path" in tmp:
                s = re.findall('d="(.*?)"', tmp)[0]
                d_path = [float(d) for d in re.findall("([-+]?\d*\.?\d*)", s) if len(d) > 0]
                oldX1, oldY1 = d_path[:2]
                a, b, c, d, e, f = re.findall('transform="matrix\((.*?),(.*?),(.*?),(.*?),(.*?),(.*?)\)', tmp)[0]

                newXs, newYs = [], []
                for chunk in chunked(d_path[2:], 6):
                    oldX = chunk[-2:][0]
                    oldY = chunk[-2:][1]
                    newX = float(a) * float(oldX) + float(c) * float(oldY) + float(e)
                    newY = float(b) * float(oldX) + float(d) * float(oldY) + float(f)
                    # if newX - newX < 1:
                    #     newX2 = newX2 - abs(newX2 - newX1) * 0.3
                    # else:
                    #     newX2 = newX2 + abs(newX2 - newX1) * 0.3
                    # if newY2 - newY1 < 1:
                    #     newY2 = newY2 - abs(newY2 - newY1) * 0.3
                    # else:
                    #     newY2 = newY2 + abs(newY2 - newY1) * 0.3
                    newXs += [newX]
                    newYs += [newY]

                newX1 = float(a) * float(oldX1) + float(c) * float(oldY1) + float(e)
                newY1 = float(b) * float(oldX1) + float(d) * float(oldY1) + float(f)

                x3 += [newX1] + newXs + [None]
                y3 += [newY1] + newYs + [None]
                print(x3)


    with open("nodes/vis/misc/tree_data_2.json", "w") as f:
        f.write(json.dumps({"x1": x1, "y1": y1, "x2": x2, "y2": y2, "x3": x3, "y3": y3}))
        f.flush()
